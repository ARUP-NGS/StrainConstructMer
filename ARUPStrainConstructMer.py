import os
from ion.plugin import IonPlugin, RunType, RunLevel, PluginCLI
import sys
import json
import simplejson
import shutil
from multiprocessing import Process, Queue
from strain_typing import Strain
from strain_typing import RunInfo
from strain_typing import SpeciesClassifier
from strain_typing import MLST
from strain_typing import AMRDatabase
from strain_typing import BaseCallerStats
from strain_typing.utility_functions import create_sql_database
from strain_typing.utility_functions import check_sql_db_for_record
from strain_typing.utility_functions import add_record_to_sql
from datetime import datetime, date

# set up for Django rendering with TSS installed apps
from django.conf import settings
from django.utils.functional import cached_property
from django.template.loader import render_to_string
from django.conf import global_settings
global_settings.LOGGING_CONFIG = None

"""
&copy; 2018, ARUP Laboratories
"""


class ARUPStrainConstructMer(IonPlugin):
    """
    This plugin performs the analysis of a sample to construct a strain suitable for comparison
    most of the strain setup should passed to the strain typing module located in this project
    This entry point handles workflow and builds the UI.
    """
    version = "1.1.0.2"  # MAJOR.MINOR.REVISION.BUILD
    DB_NAME = "STRAINS.sqlite"
    DB_DIRECTORY = "/results/plugins/scratch/"
    BACKUP_NAME = "STRAINS"
    BACKUP_DIRECTORY = None
    ENVIRONMENT = "DEV"

    runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]
    runlevels = [RunLevel.DEFAULT]
    results = ""

    @cached_property
    def startplugin_json(self):
        with open('startplugin.json', 'r') as fh:
            return json.load(fh)

    @cached_property
    def barcodes_json(self):
        with open('barcodes.json', 'r') as barcodes_handle:
            return json.load(barcodes_handle)

    # def update_progress(self, text):
    #     with open(self.startplugin['runinfo']['results_dir'] + '/progress_block.html', 'w') as html_handle:
    #         html_handle.write('<html><body>')
    #         html_handle.write("{0}".format(text))
    #         html_handle.write('</body></html>')

    def update_progress(self, text, level='info'):
        alert_hl = "alert-info"
        if level == "warn":
            alert_hl = "alert-danger"

        if text == "":
            with open(self.startplugin['runinfo']['results_dir'] + '/progress_block.html', 'w') as html_handle:
                html_handle.write("")
                return

        with open(self.startplugin['runinfo']['results_dir'] + '/progress_block.html', 'w') as html_handle:
            html_handle.write('<html><head><link rel="stylesheet" media="all" '
                              'href="/site_media/resources/bootstrap/css/bootstrap.min.css"/></head><body>')
            html_handle.write('<div class="alert {0}" role="alert">'.format(alert_hl))
            html_handle.write("<strong>{0}</strong>".format(text))
            html_handle.write('</div>')
            html_handle.write('</body></html>')

    @staticmethod
    def backup_strain(strain):
        """
        moves the plugin output directory backup directory
        :return: True if copied, False otherwise
        """
        if os.path.isdir(strain.backup_directory):
            sys.stderr.write("WARN: [{accession}, {strain_id}, {strain_instance}] "
                             "already backup up skipping\n".format(accession=strain.sample_id, strain_id=strain.sample,
                                                                   strain_instance=strain.strain_instance))
            return "ALREADY PRESENT"
        try:
            print("INFO: Coping [{0}]".format(strain.strain_directory))
            print("INFO:\tto: [{0}]".format(strain.backup_directory))
            shutil.copytree(strain.strain_directory, strain.backup_directory, symlinks=True)
            # re-adjust links to point to internal directory
            for f in os.listdir(strain.backup_directory):
                if f.startswith("IonCode_"):
                    prefix = "{acc}__{strain_id}__{instance}".format(acc=strain.sample_id, strain_id=strain.sample,
                                                                     instance=strain.strain_instance)

                    linked_name = os.path.join(strain.backup_directory, "{0}_{1}".format(prefix, f))
                    if os.path.isfile(linked_name):
                        os.unlink(linked_name)
                    os.symlink(os.path.join(strain.backup_directory, f), os.path.join(strain.backup_directory,
                                                                                      linked_name))
            return "SUCCESS"
        except OSError as e:
            raise e
        except NameError as e:
            raise e
        except:
            print("UNEXPECTED ERROR: {0}".format(sys.exc_info()[0]))
            sys.stderr.write("ERR: COULD NOT COPY TO BACKUP DRIVE\n")
            return "FAILED"

    @staticmethod
    def __get_value(key, sample_dictionary):
        value = None
        if key in sample_dictionary:
            value = sample_dictionary[key].replace(" ", "_")

        if value in [""]:
            return None
        return value

    @staticmethod
    def __is_strain_typing(value):
        if value is None:
            return False
        if value.lower() in ["st", "strtype"]:
            return True
        return False

    def validate_strains_for_analysis(self, plugin_directory, run_info=None, number_of_strains_to_process=None,
                                      basecaller_stats=None):
        """
        performs sample checks and creates list of strains to analyze

        :param plugin_directory:
        :param run_info:
        :param number_of_strains_to_process:
        :param basecaller_stats:
        :return:
        """
        sys.stderr.write("{0}\n".format("~" * 80))
        sys.stderr.write("INFO: RUNNING SAMPLE CHECKS\n")
        sys.stderr.write("{0}\n".format("~" * 80))

        strain_typing_samples = []
        for barcode in sorted(self.barcodes_json):
            metadata = self.barcodes_json[barcode]

            sample_name = self.__get_value("sample", metadata)                      # strain_id
            # sample_id = self.__get_value("sample_id", metadata)                   # accession number
            analysis_type = self.__get_value("barcode_description", metadata)       # analysis type

            is_strain_typing = self.__is_strain_typing(analysis_type)

            if sample_name is None or sample_name in ["", "None", "null"]:  # let hold onto barcode contaminates
                is_strain_typing = True
                metadata['sample'] = "{0}_{1}".format(metadata["read_count"],
                                                      metadata["bam_file"].strip("_rawlib.basecaller.bam"))
                metadata['sample_id'] = ""



            if is_strain_typing is False:  # for thermofisher deployment turn off check
                pass
            try:
                sam_id = self.startplugin_json['plan']["barcodedSamples"][metadata['sample']]["barcodeSampleInfo"][
                    barcode]['externalId'].replace(" ", "_")
            except KeyError:
                sam_id = ""
            if metadata['sample_id'] == "":
                if sam_id == "":
                    metadata['sample_id'] = ""
                else:
                    metadata['sample_id'] = sam_id
            metadata['sample'] = metadata['sample'].replace(" ", "_")
            # metadata['sample_id'] = metadata['sample_id'].replace(" ", "_").replace("-", "")

            # THESE ARE THE GUYS WE PROCESS
            if os.path.isfile(metadata['bam_filepath']):  # just check to make sure the bam is present

                already_processed, strain_instance = None, None

                # initial strain object.
                s = Strain(barcode, strain_instance, already_processed, plugin_directory, metadata=metadata,
                           run_info=run_info, startplugin_metadata=self.startplugin_json,
                           basecaller_stats=basecaller_stats)

                # create
                s.strain_directory = os.path.join(str(self.startplugin_json['runinfo']['results_dir']), barcode)
                s.add_basecalling_metadata(str(self.startplugin_json['runinfo']['basecaller_dir']))
                s.development_environment = self.ENVIRONMENT
                s.software_version = self.version

                strain_typing_samples.append(s)
                if number_of_strains_to_process is not None and len(strain_typing_samples) >= number_of_strains_to_process:  # subsampling testing
                    break

            else:
                sys.stderr.write("WARNING: Could not find {0}; skipping\n".format(metadata['bam_filepath']))

        self.update_progress("{0} samples identified for Strain Typing".format(len(strain_typing_samples)))

        sys.stderr.write("INFO: SAMPLE CHECKS DONE\n")
        sys.stderr.write("{0}\n".format("~" * 80))
        return strain_typing_samples

    def run_info_checks(self, run_info):
        """
        This method runs a set of checks on the run_info to see if the plugin should be continue
        :param run_info:
        :return: True (to launch plugin) False (stops plugin)
        """
        run_info_results = []
        sys.stderr.write("{0}\n".format("~" * 80))
        sys.stderr.write("INFO: RUNNING 'RUN INFO' TESTS\n")
        sys.stderr.write("{0}\n".format("~" * 80))

        failed_msg = "<h2>Run info Checks FAILED</h2>"

        for passed, expected_result, actual_result, test_name in run_info.check_run_info():
            run_info_results.append(passed)
            if passed is False:
                failed_msg += "<h3>Test: {0}</h3>".format(test_name)
                failed_msg += "<h3>expected result: {0}\tactual result: {1}</h3>".format(expected_result, actual_result)
                sys.stderr.write("ERROR:\tFAILED RUN INFO CHECK\n")
                sys.stderr.write("ERROR:\t{0} test FAILED\t".format(test_name))
                sys.stderr.write("[expected result: {0}||actual result: {1}\n".format(expected_result, actual_result))
            else:
                sys.stderr.write("INFO:\t{0} .... OK\t".format(test_name))
                sys.stderr.write("[expected_result: {0}||actual result: {1}]\n".format(expected_result, actual_result))

        if False in run_info_results:
            self.update_progress(failed_msg, level='warn')
            sys.stderr.write("ERROR: FINISHED RUN INFO TESTS .... ERROR\n")
            sys.stderr.write("{0}\n".format("~" * 80))
            return False
        else:
            sys.stderr.write("INFO: FINISHED RUN INFO TESTS .... OK\n")
            sys.stderr.write("{0}\n".format("~" * 80))
            return True

    def check_sql_directory(self):
        """
        check the db_location for the strain database and create it if it does not exist
        :return: the path to the strain database
        :raises: ValueError: of the db_location is not a valid path
        """
        database_file = os.path.join(self.DB_DIRECTORY, self.DB_NAME)
        if not os.path.exists(database_file):
            error_text = "INFO: Could not find the database_file at [{0}]\n".format(database_file)
            sys.stderr.write(error_text)
            self.update_progress(error_text, level='warn')

        if not os.path.isfile(database_file):
            sys.stderr.write("WARN: THE STRAIN DATABASE DOES NOT EXISTS\n")
            sys.stderr.write("INFO: CREATING DATABASE AT [{0}]\n".format(database_file))
            create_sql_database(database_file)
        else:
            sys.stderr.write("INFO: FOUND DATABASE AT [{0}]\n".format(database_file))
        return database_file

    def set_backup_and_database_directory(self):
        """
        Setup backup and database directory from global configuration

        :return: False if failed, Tru otherwise
        """
        # make sure the plugin configuration is available
        if 'pluginconfig' not in self.startplugin_json:
            print("ERROR: Failed to get plugin configuration. exiting")
            return False
        else:
            if 'DB_LOCATION' in self.startplugin_json['pluginconfig']:
                self.DB_DIRECTORY = str(self.startplugin_json['pluginconfig']['DB_LOCATION'])
                if os.path.isdir(self.DB_DIRECTORY):
                    sys.stderr.write("INFO: setting database dir to [{0}]\n".format(self.DB_DIRECTORY))
                else:
                    sys.stderr.write("ERROR: the database directory is not a valid "
                                     "directory [{0}]\n".format(self.DB_DIRECTORY))
                    self.update_progress("ERROR: the database directory is not a valid "
                                         "directory [{0}]\n".format(self.DB_DIRECTORY), level='warn')
                    return False

            # BACKUP LOCATION
            if 'STRAIN_BACKUP' in self.startplugin_json['pluginconfig']:
                if self.startplugin_json['pluginconfig']['STRAIN_BACKUP'] == 'None':
                    sys.stderr.write("INFO: No backup directory set\n")
                else:
                    self.BACKUP_DIRECTORY = self.startplugin_json['pluginconfig']['STRAIN_BACKUP']

                    if os.path.isdir(self.BACKUP_DIRECTORY):
                        sys.stderr.write("INFO: setting up backup directory to [{0}]\n".format(self.BACKUP_DIRECTORY))

                    else:
                        sys.stderr.write("ERROR: the backup directory is not a valid "
                                         "directory [{0}]\n".format(self.BACKUP_DIRECTORY))
                        self.update_progress("ERROR: the backup directory is not a valid "
                                             "directory [{0}]\n".format(self.BACKUP_DIRECTORY), level='warn')
                        return False

                    if not os.access(self.BACKUP_DIRECTORY, os.W_OK):
                        sys.stderr.write("ERROR: cannot write to the backup directory "
                                         "[{0}]\n".format(self.BACKUP_DIRECTORY))
                        self.update_progress("ERROR: cannot write to the backup directory"
                                             " [{0}]\n".format(self.BACKUP_DIRECTORY), level='warn')
                        return False

                    if not os.path.isdir(os.path.join(self.BACKUP_DIRECTORY, self.BACKUP_NAME)):
                        os.mkdir(os.path.join(self.BACKUP_DIRECTORY, self.BACKUP_NAME))
        return True

    def launch(self, data=None):
        """
        This is the main function for running samples

        :param data:
        :return:
        """
        # important directories to have handy
        results_dir = str(self.startplugin_json['runinfo']['results_dir'])  # plugin results_directory
        run_dir = str(self.startplugin_json['runinfo']['report_root_dir'])
        plugin_dir = str(self.startplugin_json['runinfo']['plugin_dir'])

        # users
        run_username = self.startplugin_json['plan']['username']
        plugin_username = self.startplugin_json['runinfo']['username']
        # self.results += str(basecaller_dir)

        basecaller_stats = BaseCallerStats(run_dir)

        # check configuration
        if not self.set_backup_and_database_directory():
            return False

        self.check_sql_directory()

        run_info = RunInfo(run_dir)
        # run_info = None


        # SET UP THE SAMPLES TO PROCESS
        strain_typing_samples = self.validate_strains_for_analysis(plugin_dir, run_info=run_info,
                                                                   number_of_strains_to_process=None,
                                                                   basecaller_stats=basecaller_stats)

        sys.stderr.write("{0}\n".format("~" * 80))
        sys.stderr.write("INFO: INITIALIZE THE RESOURCES [AMR, 16S CLASSIFIER, MLST PROFILER]\n")
        sys.stderr.write("{0}\n".format("~" * 80))

        # INITIALIZE THE AMR DATABASE
        antibiotic_gene_database = AMRDatabase(plugin_dir)

        # INITIALIZE THE CLASSIFIER DATABASE
        species_classifier = SpeciesClassifier(plugin_dir)
        # species_classifier = None

        # INITIALIZE THE MLST PROFILER
        mlst_profiler = MLST(plugin_dir)

        sys.stderr.write("{0}\n".format("~" * 80))
        sys.stderr.write("INFO: PROCESSING THE STRAINS\n")
        sys.stderr.write("{0}\n".format("~" * 80))

        # COUNT KMERS
        plugin_results = []
        for i, s in enumerate(strain_typing_samples):
            # if this returns without error a file_path will be set
            self.update_progress("Counting kmers: Strain {0} of {1}".format(i + 1, len(strain_typing_samples)))
            print("Counting kmers: Strain {0} of {1}".format(i + 1, len(strain_typing_samples)))
            s.jellyfish_file_path = s.count_kmers()

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # SUPPORT MULTIPROCESSING TO PROCESS EACH STRAIN
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cpus = 12

        if len(strain_typing_samples) < cpus:
            cpus = len(strain_typing_samples)
        self.update_progress("Setting up strains using {0} processors".format(cpus))
        q = Queue()
        jobs = strain_typing_samples
        current_processes = []
        processed_strains = []
        number_of_processes = len(jobs)
        num_of_strains_setup = 0
        for cpu in range(cpus):
            j = jobs.pop()
            p = Process(target=self.worker, args=(q, j, mlst_profiler, species_classifier, antibiotic_gene_database))
            current_processes.append(p)

        # start the current processes
        while len(current_processes) != 0:
            current_processes.pop().start()

        # keep the queue moving
        while len(jobs) != 0:
            processed_strains.append(q.get())  # pauses until job returns
            self.update_progress("INFO: {0} of {1} strains processed".format(len(processed_strains),
                                                                             number_of_processes))
            j = jobs.pop()
            num_of_strains_setup += 1
            p = Process(target=self.worker, args=(q, j, mlst_profiler, species_classifier, antibiotic_gene_database))
            p.start()

        # nothing else to start
        while number_of_processes != num_of_strains_setup:
            processed_strains.append(q.get())
            num_of_strains_setup += 1

        print("INFO: done with multiprocessing number of strains process = {0}".format(len(processed_strains)))
        q.close()
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # END OF STRAIN MULTIPROCESSING
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # SETUP OUTPUT AND UPDATE DATABASE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # set up the results so that they are sorted
        for strain in sorted(processed_strains, key=lambda x: x.name):
            sys.stderr.write("{0}\n".format("~" * 80))
            sys.stderr.write("RECORDING STRAIN INFO\n")

            # check sql to see if record exist
            already_processed, strain_instance = check_sql_db_for_record(os.path.join(self.DB_DIRECTORY, self.DB_NAME),
                                                                         strain.sample, strain.sample_id,
                                                                         strain.read_count, strain.bam_filepath)

            strain.already_processed = already_processed
            strain.strain_instance = strain_instance  # recheck version and update if needed

            strain.set_backup_directory(self.BACKUP_DIRECTORY, self.BACKUP_NAME)
            strain.create_symlinks()
            strain.write_results()

            upload_status = add_record_to_sql(os.path.join(self.DB_DIRECTORY, self.DB_NAME), strain)

            strain.upload_status = upload_status
            if self.BACKUP_DIRECTORY is None:
                sys.stderr.write("INFO: No backup made. 'Backup Directory' not set.\n")
                strain.backup_sucessful = "NOT_SET"
            else:
                strain.backup_sucessful = self.backup_strain(strain)

            print("UPLOAD_STATUS: {0}".format(upload_status))
            plugin_results.append(strain.result_summary())

        # SETUP OUTPUT PAGES
        if not settings.configured:
            settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
                               INSTALLED_APPS=('django.contrib.humanize',),
                               TEMPLATE_DIRS=(os.path.join(plugin_dir, 'templates'),))

        render_context = {
            "autorefresh": False,
            "run_name": "test_table",
            "plugin_name": os.path.basename(plugin_dir),
            "run_user_name": run_username,
            "plugin_user_name": plugin_username,
            "strain_page_link": os.path.join(self.startplugin_json['runinfo']['net_location'], 'report',
                                             str(self.startplugin_json['runinfo']['pk']), 'metal', 'plugin_out',
                                             os.path.basename(self.startplugin_json['runinfo']['plugin']
                                                              ['results_dir']), "strain_results.html"),
            "date_run": str(date.today()),
            "plugin_results": simplejson.dumps(plugin_results),
            }

        self.update_progress("")
        with open(os.path.join(results_dir, "strain_results.html"), "w") as html_fp:
            html_fp.write(render_to_string("barcode_summary.html", render_context))

        with open(os.path.join(results_dir, "status_block.html"), "w") as html_fp:
            html_fp.write(render_to_string("barcode_summary.html", render_context))
        return True

    # ##################################################################################################################
    @staticmethod
    def strain_setup(strain, mlst_profiler, species_classifier, antibiotic_gene_database):
        strain.setup_strain()
        strain.antibiotic_genes = antibiotic_gene_database.find_antibiotic_genes(strain)
        strain.write_antibiotic_genes()
        strain.mlst_profiles, profiles_loaded = mlst_profiler.find_mlst_profiles(strain)
        strain.hits = species_classifier.compare_references(strain.kmers)
        strain.write_hits()
        strain.write_mlst()
        strain.write_hits_to_html()
        strain.create_resistant_table()
        return strain

    # THIS IS THE METHOD THAT DOES MULTIPROCESSING
    def worker(self, queue, strain, mlst_profiler, species_classifier, antibiotic_gene_database):
        queue.put(self.strain_setup(strain, mlst_profiler, species_classifier, antibiotic_gene_database))


if __name__ == "__main__":
    PluginCLI()
