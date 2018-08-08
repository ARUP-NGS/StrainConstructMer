import os
import sys
import json
import gzip
from subprocess import Popen
from subprocess import PIPE
from subprocess import check_call
from subprocess import CalledProcessError
from collections import OrderedDict
from .plots import produce_histograms
from .utility_functions import reverse_complement
from .utility_functions import timeit
import datetime
import simplejson
from Bio import SeqIO
from django.conf import settings
from django.template.loader import render_to_string
from django.conf import global_settings
global_settings.LOGGING_CONFIG = None


class Strain(object):
    """
    This is the primary object of the project and represents.

    """

    jellyfish_app = "/results/plugins/scratch/apps/jellyfish-2.2.7/bin/jellyfish"
    KMER_SIZE = 31

    # QUALITY CUTOFFS
    COVERAGE_CUTOFF = 25
    GENOME_SIZE_CUTOFF = 1000000
    DIFFERENCE_FROM_EXPECTED_CUTOFF = 15
    PCT_QC20_BASES_CUTOFF = 70
    COVERAGE_CAP_READS = 10000000

    def __init__(self, name, strain_instance, already_processed, plugin_directory, metadata=None,
                 startplugin_metadata=None, run_info=None, is_reference=False, basecaller_stats=None):

        self.basecaller_stats = basecaller_stats
        self.name = name
        self.strain_instance = strain_instance
        self.software_version = None
        self.already_processed = already_processed
        self.plugin_dir = plugin_directory
        self.is_reference = is_reference
        self.metadata = metadata
        self.start_plugin_metadata = startplugin_metadata

        # properties
        self.backup_sucessful = False
        self.control_sample = False
        self.backup_directory = None
        self.upload_status = None
        self.development_environment = None
        self.server = None
        self.date_created = datetime.datetime.today().isoformat()
        self._read_count = None
        self._barcode_adapter = None
        self._sample = None
        self._accession = None
        self._strain_id = None
        self._strain_directory = None
        self._bam_filepath = None
        self._bam_file = None
        self._barcode_name = None
        self._barcode_sequence = None
        self._barcode_index = None
        self._sample_id = None
        self.report_directory = None
        self._rrna_reads_path = None
        self._rrna_contig_path = None
        self._full_length_ssu_path = None
        self.antibiotic_genes = None
        self.mlst_profile_path = None
        self.named_strain_directory = None

        # web icons
        self.failed_icon = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/cancel-square.png"
        self.passed_icon = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/accept-square.png"
        self.alert_icon = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/alert-triangle-yellow.png"
        self.orange_dot = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/warning.png"
        self.green_dot = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/success_checkmark.png"
        self.yellow_dot = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/warning.png"
        self.red_dot = self.plugin_dir.replace("/results", "") + "/pluginMedia/images/failed.png"

        # basecaller info
        self.bc_max_read_length_q17 = None
        self.bc_mean_read_length_q17 = None
        self.bc_read_count_q17 = None
        self.bc_base_count_q17 = None
        self.bc_max_read_length_q20 = None
        self.bc_mean_read_length_q20 = None
        self.bc_read_count_q20 = None
        self.bc_base_count_q20 = None
        self.bc_max_read_length = None
        self.bc_mean_read_length = None
        self.bc_read_count = None
        self.bc_base_count = None
        self.histo_plot_path = None
        self.reads_to_analyze = None
        self.kmer_histo = None
        self.antibiotic_path = None
        self.strain_json_path = None

        # run_info object
        self.run_info = run_info

        # measured attributes
        self._jellyfish_file_path = None
        self._jellyfish_histo_path = None
        self._histo_valley_index = None
        self._histo_valley_freq = None
        self._histo_valley_count = None
        self._kmer_count = None
        self._distinct_kmers = None
        self._coverage = None
        self.classifier_html = None
        self.classifier_text_path = None

        # kmers
        self._kmer_path = None
        self._kmers = None
        self._hits = None
        self._species_predicted_genome_size = None
        self._genus_predicted_genome_size = None

        self.mlst_profiles = None

        self.passed = False

        # CHEF INFO
        self.chef_info = {
            "Ion S5 Cleaning Solution": {
                "partNumber": None,
                "lotNumber": None,
                "expDate": None
            },
            "Ion S5 Sequencing Reagents": {
                "partNumber": None,
                "lotNumber": None,
                "expDate": None
            },
            "Ion S5 Wash Solution": {
                "partNumber": None,
                "lotNumber": None,
                "expDate": None
            }
        }

        self.__parse_metadata()

    # def parse_initLog(self, report_directory):
    #
    #     choice = ["Ion S5 Cleaning Solution",
    #               "Ion S5 Wash Solution",
    #               "Ion S5 Sequencing Reagents"]
    #
    #     if self.run_info is None:
    #         return None
    #
    #     init_path = os.path.join(report_directory, "InitLog.txt")
    #     if not os.path.isfile(init_path):
    #         init_path = os.path.join(self.run_info.experiment_directory, "InitLog.txt")
    #
    #     if not os.path.isfile(init_path):
    #         return None
    #
    #     current_choice = None
    #     for line in open(init_path):
    #         line = line.strip().split(":")
    #         if len(line) != 2:
    #             continue
    #         key, value = line[0].strip(), line[1].strip()
    #         if key == "productDesc":
    #             if value in choice:
    #                 current_choice = value
    #         if current_choice is None:
    #             continue
    #         if key in ["partNumber", "lotNumber", "expDate"]:
    #             self.chef_info[current_choice][key] = value
    #     return None

    def __create_symlink_strain_directory(self):
        """
        this creates the symlink for the directory

        :return:
        """
        clean_sample_name = self.sample.replace(" ", "_")
        if os.path.isdir(self.strain_directory):
            named_strain_directory = "{0}_{1}".format(clean_sample_name, os.path.basename(self.strain_directory))
            os.symlink(self.strain_directory, os.path.join(os.path.dirname(self.strain_directory),
                                                           named_strain_directory))
            self.named_strain_directory = os.path.join(os.path.dirname(self.strain_directory),
                                                       named_strain_directory)

    def create_symlinks(self):
        """
        function creates symlinks with the sample name as the prefix
        if any of the directory names are none it will skip


        :return: True if symlinks create, False otherwise
        """
        # only perform if the sample name is valid
        if self.sample is None or self.sample == "":
            sys.stderr.write("no symlinks created sample name empty of None\n")
            return False

        # create a symlink directory
        self.__create_symlink_strain_directory()
        # clean_sample_name = self.sample.replace(" ", "_")

        if self.sample_id in [None, ""]:
            prefix = "{sample}__{version}".format(sample=self.sample, version=self.strain_instance)
        else:
            prefix = "{acc}__{strain_id}__{version}".format(acc=self.sample_id, strain_id=self.sample,
                                                            version=self.strain_instance)

        # try to symlink these files
        paths_to_link = [self.rrna_reads_path,
                         self.rrna_contig_path,
                         self.full_length_ssu_path,
                         self.mlst_profile_path,
                         self.kmer_path,
                         self.jellyfish_file_path,
                         self.jellyfish_histo_path,
                         self.classifier_html,
                         self.classifier_text_path,
                         self.antibiotic_path,
                         self.histo_plot_path,
                         self.strain_json_path]

        # link all of our data
        for p in paths_to_link:
            if p is not None and os.path.isfile(p):
                os.symlink(p, os.path.join(os.path.dirname(p), "{0}_{1}".format(prefix,
                                                                                os.path.basename(p))))

        # link the bam file
        bl = os.path.join(self.strain_directory, self.bam_file)
        if os.path.isfile(self.bam_filepath):
            os.symlink(self.bam_filepath, bl)
            os.symlink(bl, os.path.join(os.path.dirname(bl), "{0}_{1}".format(prefix, os.path.basename(bl))))

    @timeit
    def sequester_rrna_reads(self, read_cutoff=500, kmer_cutoff=20, bam_file=True, gzipped=False):  # pragma: no cover
        """
        This function creates a system call to the sequester app and creates a results file
        samtools bam2fq IonCode_0117_rawlib.basecaller.bam | head | ./query_per_sequence SSU.jf /dev/fd/0
        :return: path_to_rrna_reads
        """
        rrna_reads = os.path.join(self.strain_directory, self.name + "_16S_reads.fq")

        with open(rrna_reads, "w") as wf:
            if bam_file:
                bam_to_fastq_process = Popen(["bamToFastq", "-i", self.bam_filepath, "-fq", "/dev/stdout"], stdout=PIPE)
                sequester_process = Popen([self.sequester_app, self.sequester_16S_db, str(read_cutoff),
                                           str(kmer_cutoff), "1", "/dev/fd/0"], stdin=bam_to_fastq_process.stdout,
                                          stdout=PIPE, stderr=PIPE)
                bam_to_fastq_process.stdout.close()
                for line in sequester_process.stdout:
                    wf.write(line)

                for line in sequester_process.stderr:
                    sys.stderr.write("INFO: [sequester] {0}".format(line))

                sequester_process.wait()

            else:
                if gzipped:
                    gzip_process = Popen(["gzip", "--decompress", "--to-stdout", self.bam_filepath], stdout=PIPE)

                    sequester_process = Popen([self.sequester_app, self.sequester_16S_db, str(read_cutoff),
                                              str(kmer_cutoff), "1", "/dev/fd/0"], stdin=gzip_process.stdout,
                                              stdout=PIPE, stderr=PIPE)

                    for line in sequester_process.stdout:
                        wf.write(line)

                    for line in sequester_process.stderr:
                        sys.stderr.write("INFO: [sequester] {0}".format(line))

                    sequester_process.wait()
                    sequester_process.stdout.close()
                    gzip_process.stdout.close()

                else:
                    sequester_process = Popen([self.sequester_app, self.sequester_16S_db, str(read_cutoff),
                                               str(kmer_cutoff), "1", self.bam_filepath], stdout=PIPE)
                    for line in sequester_process.stdout:
                        wf.write(line)
                    sequester_process.wait()

        self.rrna_reads_path = rrna_reads

    def write_mlst(self):  # pragma: no cover
        """
        writes mlst results to strain_directory

        :return: None
        """
        self.mlst_profile_path = os.path.join(self.strain_directory, self.name + "_mlst_profiles.txt")
        with open(self.mlst_profile_path, "w") as wf:
            if self.mlst_profiles is None:
                return
            for i, mlst in enumerate(sorted(self.mlst_profiles, key=lambda x: x.average_coverage, reverse=True)):
                wf.write("{0}".format(mlst.print_profile()))
        return

    def write_hits_to_html(self):  # pragma: no cover
        """
        writes an HTML file to strain_directory containing the classifier hits
        :return: None
        """
        if not settings.configured:
            settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
                               INSTALLED_APPS=('django.contrib.humanize',),
                               TEMPLATE_DIRS=(os.path.join(self.plugin_dir, 'templates'),))

        formatted_hits = []
        if self.hits is not None:
            for h in self.hits:
                _d = {"predicted_size": self.distinct_kmers}
                for k, v in h.__dict__.items():
                    if k == "kmer_identity":
                        v = float("{0:.2f}".format(float(v)))
                        _d[k] = v
                    else:
                        _d[k] = v
                formatted_hits.append(_d)

        if self.start_plugin_metadata is None:
            run_user_name = "CLI"
            plugin_user_name = "CLI"
            pk = ""
            self.classifier_html = os.path.join(self.strain_directory, self.name + "_classifier_hits.html")
        else:
            run_user_name = self.start_plugin_metadata['plan']['username']
            plugin_user_name = self.start_plugin_metadata['runinfo']['username']
            pk = "/" + os.path.join("report", str(self.start_plugin_metadata['runinfo']['pk'])) + "/#plugins-section"
            server = self.start_plugin_metadata['runinfo']['net_location']
            report_number = str(self.start_plugin_metadata['runinfo']['pk'])
            results_dir = os.path.basename(self.start_plugin_metadata['runinfo']['plugin']['results_dir'])
            self.classifier_html = os.path.join(server, 'report', report_number, 'metal', 'plugin_out', results_dir,
                                                self.name, self.name + "_classifier_hits.html")
        render_context = {
            "autorefresh": False,
            "run_name": "test_table",
            "run_user_name": run_user_name,
            "plugin_user_name": plugin_user_name,
            "date_run": str(datetime.date.today()),
            "pk": pk,
            "barcode": self.name,
            "sample": self.sample,
            "strain_results": simplejson.dumps(formatted_hits),
            # "help_dict": help_dictionary,
        }

        with open(os.path.join(self.strain_directory, self.name + "_classifier_hits.html"), "w") as html_fp:
            html_fp.write(render_to_string("classifier_table.html", render_context))

    def write_hits(self):  # pragma: no cover
        """
        Writes a text file to strain directory containing the classifier hits
        :return: None
        """
        self.classifier_text_path = os.path.join(self.strain_directory, self.name + "_classifier_hits.txt")
        with open(self.classifier_text_path, "w") as wf:
            if self.hits is None:
                return
            for i, hit in enumerate(self.hits):
                header, row = hit.print_hit(delimiter=",")
                if i == 0:
                    wf.write("{0}\n".format(header))
                wf.write("{0}\n".format(row))
        return

    def write_antibiotic_genes(self):  # pragma: no cover
        """
        writes antibiotic resistant hits to json formatted file in strain_directory
        :return: None
        """
        self.antibiotic_path = os.path.join(self.strain_directory, self.name + "_full_antibiotic_genes.json")
        with open(self.antibiotic_path, "w") as wf:
            if self.antibiotic_genes is None:
                return
            json.dump(self.antibiotic_genes, wf, indent=4)

    @staticmethod
    def __split_species(species_name):
        """
        helper method that splits the species name and returns a short strain_instance
        (e.g. Staphylococcus aureus; S. aureus)

        :param species_name: Str: species name
        :return: Str: abbv species name
        """
        genus, species = species_name.split(" ")
        return "{0}. {1}".format(genus[0], species)

    def __passed_coverage(self, coverage_cutoff=None):
        if coverage_cutoff is None:
            coverage_cutoff = self.COVERAGE_CUTOFF

        if self.coverage < coverage_cutoff:
            return False
        else:
            return True

    def __passed_q20(self, pct_qc_bases=None):
        if pct_qc_bases is None:
            pct_qc_bases = self.PCT_QC20_BASES_CUTOFF

        qc_bases = self.__get_pct_q20_bases()
        if qc_bases is not None and qc_bases < pct_qc_bases:
            return False
        else:
            return True

    def __passed_genome_size(self, genome_size_cutoff=None):
        if genome_size_cutoff is None:
            genome_size_cutoff = self.GENOME_SIZE_CUTOFF

        if self.distinct_kmers < genome_size_cutoff:
            return False
        else:
            return True

    def get_genome_min_genome_size(self):
        if self.hits is None or len(self.hits) == 0:
            return None
        if self.hits[0].min_genome_size_species is not None:
            return self.hits[0].min_genome_size_species
        if self.hits[0].min_genome_size_genus is not None:
            return self.hits[0].min_genome_size_genus
        return None

    def get_genome_max_genome_size(self):
        if self.hits is None or len(self.hits) == 0:
            return None
        if self.hits[0].max_genome_size_species is not None:
            return self.hits[0].max_genome_size_species
        if self.hits[0].max_genome_size_genus is not None:
            return self.hits[0].max_genome_size_genus
        return None

    def __passed_genome_range(self):
        if self.hits is None or len(self.hits) == 0:
            return None
        if self.distinct_kmers < (self.get_genome_min_genome_size() * 0.9):
            return False
        if self.distinct_kmers > (self.get_genome_max_genome_size() * 1.1):
            return False

        return True

    def strain_passed_qc(self, coverage_cutoff=None, genome_size_cutoff=None, pct_qc_bases=None):
        """
        returns true if the strain passes set QC parameters

        :param coverage_cutoff: The minimum coverage allowed
        :param genome_size_cutoff: The minimum genomes size allowed
        :param pct_qc_bases: The percentage of q20 bases acceptable
        :return: Bool True if passes QC
        """
        passed_coverage = self.__passed_coverage(coverage_cutoff=coverage_cutoff)
        passed_q20 = self.__passed_q20(pct_qc_bases=pct_qc_bases)
        passed_genome_size = self.__passed_genome_size(genome_size_cutoff=genome_size_cutoff)
        passed_genome_range = self.__passed_genome_range()

        if False in [passed_coverage, passed_q20, passed_genome_size, passed_genome_range]:
            return False
        else:
            return True

    def __get_report_directory(self):  # pragma: no cover
        """
        Try to return the report directory url based on plugin metadata data

        # http://s5xldev/report/4/metal/plugin_out/PluginStrainTyping_out.264/IonCode_0117
        :return: None if no plugin data or a str representing url
        """
        if self.start_plugin_metadata is None:
            return None
        server = self.start_plugin_metadata['runinfo']['net_location']
        report_number = str(self.start_plugin_metadata['runinfo']['pk'])
        results_dir = os.path.basename(self.start_plugin_metadata['runinfo']['plugin']['results_dir'])
        report_directory = os.path.join(server, 'report', report_number, 'metal', 'plugin_out', results_dir, self.name)
        return report_directory

    def set_backup_directory(self, base_location, dir_name):

        if base_location is None:
            return None

        run_date = datetime.datetime.strptime(self.run_info.date, "%Y-%m-%dT%H:%M:%S")
        year_folder = str(run_date.year)

        month_folder = "{0:02}_{1}".format(run_date.month, run_date.strftime("%B"))

        if not os.path.exists(os.path.join(base_location, dir_name, year_folder)):
            os.mkdir(os.path.join(base_location, dir_name, year_folder))

        if not os.path.exists(os.path.join(base_location, dir_name, year_folder, month_folder)):
            os.mkdir(os.path.join(base_location, dir_name, year_folder, month_folder))

        if self.sample_id in [None, ""]:
            prefix = "{sample}__{version}".format(sample=self.sample, version=self.strain_instance)
        else:
            prefix = "{acc}__{strain_id}__{version}".format(acc=self.sample_id, strain_id=self.sample,
                                                            version=self.strain_instance)
        self.backup_directory = os.path.join(base_location, dir_name, year_folder, month_folder, prefix)

    def prepare_metadata_output(self):
        if self.start_plugin_metadata is None:
            return ""
        _d = self.start_plugin_metadata
        _d.pop("plan")
        return _d

    def write_results(self):
        """
        writes a json formatted file to the strain_directory.
        :return: the write dictionary containing output values
        """
        write_dictionary = OrderedDict()
        write_dictionary['accession'] = self.accession
        write_dictionary['strain_id'] = self.strain_id
        write_dictionary["strain_instance"] = self.strain_instance
        write_dictionary["software_version"] = self.software_version
        write_dictionary["control_sample"] = self.control_sample
        write_dictionary["strain_directory"] = self.strain_directory
        write_dictionary["backup_directory"] = self.backup_directory
        # write_dictionary["backup_sucessful"] = self.backup_sucessful

        write_dictionary["development_environment"] = self.development_environment

        if self.start_plugin_metadata is not None:
            write_dictionary["server"] = self.start_plugin_metadata['runinfo']['net_location']
            write_dictionary["run_user_name"] = self.start_plugin_metadata['plan']['username']
            write_dictionary["plugin_user_name"] = self.start_plugin_metadata['runinfo']['username']

        write_dictionary["barcode"] = self.name
        write_dictionary["sample"] = self.sample
        write_dictionary["sample_id"] = self.sample_id

        write_dictionary["processed_before"] = self.already_processed

        write_dictionary["date_created"] = self.date_created
        write_dictionary["plugin_directory"] = self.plugin_dir

        # write_dictionary["is_reference"] = self.is_reference
        write_dictionary["metadata"] = self.metadata
        write_dictionary["read_count"] = self.read_count
        write_dictionary["barcode_adapter"] = self.barcode_adapter

        write_dictionary["bam_filepath"] = self._bam_filepath
        write_dictionary["bam_file"] = self._bam_file
        write_dictionary["barcode_name"] = self.barcode_name
        write_dictionary["barcode_sequence"] = self.barcode_sequence
        write_dictionary["barcode_index"] = self.barcode_index

        write_dictionary["plugin_metadata"] = self.prepare_metadata_output()
        write_dictionary["report_directory"] = self.report_directory
        write_dictionary["count_intersection"] = self.histo_valley_count
        write_dictionary["frequency_cutoff"] = self.histo_valley_freq  # same as above

        write_dictionary["max_read_length_q20"] = self.bc_max_read_length_q20
        write_dictionary["mean_read_length_q20"] = self.bc_mean_read_length_q20
        write_dictionary["read_count_q20"] = self.bc_read_count_q20
        write_dictionary["base_count_q20"] = self.bc_base_count_q20

        write_dictionary["max_read_length"] = self.bc_max_read_length
        write_dictionary["mean_read_length"] = self.bc_mean_read_length
        write_dictionary["read_count"] = self.bc_read_count
        write_dictionary["base_count"] = self.bc_base_count
        write_dictionary["pct_q20_bases"] = self.__get_pct_q20_bases()
        write_dictionary["reads_analyzed"] = self.reads_to_analyze
        write_dictionary["sample"] = self.sample
        write_dictionary["coverage"] = self.coverage
        write_dictionary["genome_size"] = self.distinct_kmers
        write_dictionary["difference_from_expected"] = self.__calculate_genome_size_difference()

        write_dictionary["histo_image"] = self.histo_plot_path
        write_dictionary["passed_qc"] = self.strain_passed_qc()
        write_dictionary["passed_coverage_cutoff"] = self.__passed_coverage()
        write_dictionary["coverage_cutoff"] = self.COVERAGE_CUTOFF
        write_dictionary["passed_genome_size_cutoff"] = self.__passed_genome_size()
        write_dictionary["genome_size_cutoff"] = self.GENOME_SIZE_CUTOFF
        write_dictionary["passed_q20_cutoff"] = self.__passed_q20()
        write_dictionary["q20_cutoff"] = self.PCT_QC20_BASES_CUTOFF
        write_dictionary["passed_genome_range"] = self.__passed_genome_range()
        write_dictionary["min_known_genome_size"] = self.get_genome_min_genome_size()
        write_dictionary["max_known_genome_size"] = self.get_genome_max_genome_size()
        write_dictionary["chef_info"] = self.chef_info

        if self.run_info is None:
            write_dictionary["run_info"] = None
        else:
            write_dictionary["run_info"] = self.run_info.__dict__

        # format classifier hits
        write_dictionary['classifier_hits'] = []
        if self.hits is not None:
            for h in self.hits:
                write_dictionary['classifier_hits'].append(h.formatted_output())

        # format mlst hits
        write_dictionary["mlst_profiles"] = []
        if self.mlst_profiles is not None:
            for profile in sorted(self.mlst_profiles, key=lambda x: x.average_coverage, reverse=True):
                write_dictionary["mlst_profiles"].append(profile.__dict__)

        # format antibiotic hits:
        write_dictionary["antibiotic_genes"] = []
        gene_list = set([])
        if self.antibiotic_genes is not None:
            for a_gene in self.antibiotic_genes:
                if a_gene['gene'] not in gene_list:
                    a_gene.pop("coverage")
                    write_dictionary["antibiotic_genes"].append(a_gene)
                    gene_list.add(a_gene['gene'])

        write_dictionary["kmer_path"] = self.kmer_path
        write_dictionary['species_predicted_genome_size'] = self.species_predicted_genome_size
        write_dictionary['genus_predicted_genome_size'] = self.genus_predicted_genome_size

        self.strain_json_path = os.path.join(self.strain_directory, self.name + "_strain.json")
        with open(self.strain_json_path, "w") as json_output:
            json.dump(write_dictionary, json_output, indent=4)
        return write_dictionary

    def result_summary(self):
        """
        creates a dictionary to that is passed to the html renderer
        return dictionary
        """
        results_dictionary = dict()
        results_dictionary["plugin_directory"] = self.plugin_dir
        results_dictionary["is_reference"] = self.is_reference
        results_dictionary["metadata"] = self.metadata
        results_dictionary["strain_instance"] = self.strain_instance
        results_dictionary["read_count"] = self.read_count
        results_dictionary["barcode_adapter"] = self.barcode_adapter
        results_dictionary["strain_directory"] = self.strain_directory
        results_dictionary["bam_filepath"] = self._bam_filepath
        results_dictionary["bam_file"] = self._bam_file
        results_dictionary["barcode_name"] = self.barcode_name
        results_dictionary["barcode_sequence"] = self.barcode_sequence
        results_dictionary["classifier_html"] = self.classifier_html
        results_dictionary["barcode_index"] = self.barcode_index
        results_dictionary["accession"] = self.accession
        results_dictionary["strain_id"] = self.strain_id
        results_dictionary["sample_id"] = self.sample_id
        results_dictionary["kmer_path"] = self.kmer_path
        results_dictionary["report_directory"] = self.report_directory
        results_dictionary["count_intersection"] = self.histo_valley_count
        results_dictionary["frequency_cutoff"] = self.histo_valley_freq  # same as above

        # basecaller info
        results_dictionary["max_read_length_q17"] = self.bc_max_read_length_q17
        results_dictionary["mean_read_length_q17"] = self.bc_mean_read_length_q17
        results_dictionary["read_count_q17"] = self.bc_read_count_q17
        results_dictionary["base_count_q17"] = self.bc_base_count_q17

        results_dictionary["max_read_length_q20"] = self.bc_max_read_length_q20
        results_dictionary["mean_read_length_q20"] = self.bc_mean_read_length_q20
        results_dictionary["pct_q20_bases"] = self.__get_pct_q20_bases()

        results_dictionary["read_count_q20"] = self.bc_read_count_q20
        results_dictionary["base_count_q20"] = self.bc_base_count_q20

        results_dictionary["max_read_length "] = self.bc_max_read_length
        results_dictionary["mean_read_length "] = self.bc_mean_read_length
        results_dictionary["read_count "] = self.bc_read_count
        results_dictionary["base_count "] = self.bc_base_count
        results_dictionary["reads_analyzed"] = self.reads_to_analyze
        results_dictionary["barcode"] = self.name.split("_")[-1]
        results_dictionary["coverage_format"] = "Green" if self.__passed_coverage() else "Red"
        results_dictionary["genome_size_format"] = "Green" if self.__passed_genome_size() else "Red"
        results_dictionary["genome_size_icon"], results_dictionary["genome_size_help"] = self.get_genome_size_icon()
        results_dictionary["mongo_icon"], results_dictionary["mongo_help"] = self.get_mongo_icon()
        results_dictionary["backup_icon"], results_dictionary["backup_help"] = self.get_backup_icon()

        results_dictionary["qc_base_format"] = "Green" if self.__passed_q20() else "Red"
        results_dictionary["sample"] = self.sample
        results_dictionary["coverage"] = self.coverage
        results_dictionary["genome_size"] = self.distinct_kmers
        results_dictionary["difference_from_expected"] = self.__calculate_genome_size_difference()
        results_dictionary["chef_info"] = self.chef_info

        if self.histo_plot_path is None:
            results_dictionary["histo_image"] = None
        else:
            results_dictionary["histo_image"] = self.histo_plot_path.replace("/results/analysis", "")
        results_dictionary["passed_qc"] = self.strain_passed_qc()

        results_dictionary["qc_icon"] = self.get_qc_icon()

        results_dictionary["count_intersection"] = self.histo_valley_count
        if self.hits is None:
            results_dictionary["top_hit_kid"] = "NONE"
            results_dictionary["top_hit_species"] = "NONE"
            results_dictionary["top_hit_accession"] = "NONE"
        else:
            try:
                results_dictionary["top_hit_kid"] = float("{0:.1f}".format(self.hits[0].kmer_identity))
                results_dictionary["top_hit_species"] = self.__split_species(self.hits[0].species)
                results_dictionary["top_hit_accession"] = self.hits[0].accession.split(".")[0]
            except IndexError:
                results_dictionary["top_hit_kid"] = "NONE"
                results_dictionary["top_hit_species"] = "NONE"
                results_dictionary["top_hit_accession"] = "NONE"

        if self.mlst_profiles is None:
            results_dictionary["mlst_profile"] = "NONE"
        else:
            s = ""
            for p in sorted(self.mlst_profiles, key=lambda x: x.average_coverage, reverse=True):
                s += "{0}\n".format(str(p).strip())
            results_dictionary["mlst_profile"] = s
        return results_dictionary

    def get_genome_size_icon(self):
        if self.__calculate_genome_size_difference() is None:
            return self.orange_dot, "The expected genome size is unknown"

        elif self.__passed_genome_range():
            return (self.green_dot, "The expected genome size is within 10% of the min and max genome size "
                                    "observed for this species/genus")
        else:
            return (self.red_dot, "The expected genome size is NOT within 10% of the min and max genome size "
                                  "observed for this species/genus")

    def get_mongo_icon(self):
        if self.upload_status is None or self.already_processed:
            return self.green_dot, "The sample is already in the database; replaced old record"

        if self.upload_status is "FAILED":
            return self.red_dot, "An error occured when attempting to write to database"

        if self.upload_status is "SUCCESS":
            return self.green_dot, "The sample was successfully backup"

        return self.red_dot, "FAILED: unknown error"

    def get_backup_icon(self):
        if self.backup_sucessful == "ALREADY PRESENT":
            return self.orange_dot, "The sample has been backed up previously"
        if self.backup_sucessful == "FAILED":
            return self.red_dot, "An error occurred backing up the sample"
        if self.backup_sucessful == "SUCCESS":
            return self.green_dot, "Strain backup to: {0}".format(self.backup_directory)
        if self.backup_sucessful == "NOT_SET":
            return self.orange_dot, "Backup directory not set."
        return self.orange_dot, "unknown error"

    def get_qc_icon(self):
        """
        helper method to get the correct result icon
        :return: path to icon
        """
        if self.accession == "BARCODE_CONTAMINATE":
            return self.alert_icon
        if self.strain_passed_qc():
            return self.passed_icon
        else:
            return self.failed_icon

    def __get_pct_q20_bases(self):
        """
        helper method to get the pct of q20 basses
        :return:
        """
        if self.bc_base_count_q20 is None or self.bc_base_count is None:
            return None
        return float("{0:.1f}".format(float(self.bc_base_count_q20) / float(self.bc_base_count) * 100))

    def __repr__(self):  # pragma: no cover
        """
        represtation of the strain object
        :return: string
        """
        s = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        s += "name:\t{0}\n".format(self.name)
        s += "read_count:\t{0}\n".format(self.read_count)
        s += "barcode_adapter:\t{0}\n".format(self.barcode_adapter)
        s += "sample:\t{0}\n".format(self.sample)
        s += "strain_dir:\t{0}\n".format(self.strain_directory)
        s += "bam_path:\t{0}\n".format(self.bam_filepath)
        s += "bam_file:\t{0}\n".format(self.bam_file)
        s += "barcode_name:\t{0}\n".format(self.barcode_name)
        s += "barcode_sequence:\t{0}\n".format(self.barcode_sequence)
        s += "barcode_index:\t{0}\n".format(self.barcode_index)
        s += "sample_id:\t{0}\n".format(self.sample_id)
        s += "mean_read_length [all, q17, q20]:\t[{0}, {1}, {2}]\n".format(self.bc_mean_read_length,
                                                                           self.bc_mean_read_length_q17,
                                                                           self.bc_mean_read_length_q20)
        s += "base_count [all, q17, q20]:\t[{0}, {1}, {2}]\n".format(self.bc_base_count,
                                                                     self.bc_base_count_q17,
                                                                     self.bc_base_count_q20)
        s += "read_count [all, q17, q20]:\t[{0}, {1}, {2}]\n".format(self.bc_read_count,
                                                                     self.bc_read_count_q17,
                                                                     self.bc_read_count_q20)

        if self.reads_to_analyze is not None:
            s += "reads_analyzed:\t{0:,}\n".format(min(self.reads_to_analyze, self.read_count))

        s += "jellyfish_path:\t{0}\n".format(self.jellyfish_file_path)
        s += "jellyfish_histo_path:\t{0}\n".format(self.jellyfish_histo_path)
        s += "histo_vally [index, freq, count]:\t[{0}, {1}, {2}]\n".format(self.histo_valley_index,
                                                                           self.histo_valley_freq,
                                                                           self.histo_valley_count)
        s += "total observed kmer [above cutoff]:\t{0}\n".format(self.kmer_count)
        s += "distinct_kmers (e.g. genome_size):\t{0}\n".format(self.distinct_kmers)
        s += "coverage:\t{0}\n".format(self.coverage)
        if self.__calculate_genome_size_difference() is not None:
            s += "difference from expected genome:\t{0}\n".format(self.__calculate_genome_size_difference())
        else:
            s += "difference from expected genome:\tNA\n"

        s += "species hits:\n"
        if self.hits is not None:
            for h in self.hits:
                s += "\t{0}\n".format(str(h))
        else:
            s += "\tNone\n"

        s += "mlst profiles:\n"
        if self.mlst_profiles is not None:
            for p in self.mlst_profiles:
                s += "\t{0}\n".format(p)
        else:
            s += "\tNone\n"

        s += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return s

    def __calculate_genome_size_difference(self):
        """
        returns the difference of genomes size of the top hit

        :return: float with one decimal precision if value can be calulated None otherwise
        """
        if self.distinct_kmers is None:
            return None

        elif self.species_predicted_genome_size is not None:
            v = float(abs(
                self.species_predicted_genome_size - self.distinct_kmers)) / self.species_predicted_genome_size * 100.0
            return float("{0:.1f}".format(v))
        elif self.genus_predicted_genome_size is not None:
            v = float(abs(
                self.genus_predicted_genome_size - self.distinct_kmers)) / self.genus_predicted_genome_size * 100.0
            return float("{0:.1f}".format(v))
        else:
            return None

    @timeit
    def setup_strain(self):
        """
        Method can be run after the kmers have been counted to generate strain calculations
        :return: None
        """
        self.__generate_kmer_histo()
        self.histo_valley_index, self.histo_valley_freq, self.histo_valley_count = self.__get_initial_shallow()
        self.kmer_path, self.kmers = self.__load_kmers()
        total_kmers = 0
        distinct_kmers = 0
        for freq, count in self.kmer_histo[self.histo_valley_index:]:
            total_kmers += freq * count
            distinct_kmers += count

        self.distinct_kmers = distinct_kmers
        self.kmer_count = total_kmers
        if total_kmers == 0 or distinct_kmers == 0:
            self.coverage = 0
            histo_path = produce_histograms([(0, 0,), (1, 1)], self.coverage, self.sample, self.histo_valley_freq,
                                            self.histo_valley_count, self.name, self.distinct_kmers,
                                            os.path.join(self.strain_directory, self.name + "_histo.png"))

            if histo_path is None:
                self.histo_plot_path = self.plugin_dir.replace("/results", "") + \
                                       "/pluginMedia/images/cancel-square-64.png"
            else:
                self.histo_plot_path = histo_path

        else:
            self.coverage = int(float(total_kmers) / float(distinct_kmers))

            histo_path = produce_histograms(self.kmer_histo, self.coverage, self.sample, self.histo_valley_freq,
                                            self.histo_valley_count, self.name, self.distinct_kmers,
                                            os.path.join(self.strain_directory, self.name + "_histo.png"))

            if histo_path is None:
                self.histo_plot_path = self.plugin_dir.replace("/results", "") + \
                                       "/pluginMedia/images/cancel-square-64.png"
            else:
                self.histo_plot_path = histo_path

        self.report_directory = self.__get_report_directory()

    def __load_kmers(self):
        """
        helper function that loads the kmers from jellyfish output

        :return: kmer_output_path, dictionary of kmer: frequency
        """
        if os.path.isfile(self.jellyfish_file_path) is False:
            raise ValueError("Could not find jellyfish file to generate kmers "
                             "[{0}]".format(self.jellyfish_file_path))
        kmer_outfile = os.path.join(self.strain_directory, self.name + ".kmers")

        frequency_cutoff = 3
        if self.histo_valley_freq - 5 > frequency_cutoff:
            frequency_cutoff = self.histo_valley_freq - 5

        try:
            check_call([self.jellyfish_app, "dump", "--output", kmer_outfile, "--column", "--tab",
                        "--lower", str(frequency_cutoff), self.jellyfish_file_path])
            check_call(["gzip", kmer_outfile])

            kmer_outfile = kmer_outfile + ".gz"

        except CalledProcessError as err:
            print(err)

        if os.path.isfile(kmer_outfile) is False:  # there were no kmers above cutoff
            self.kmers = {}
            return None

        kmer_dictionary = {}  # passing kmers
        for line in gzip.open(kmer_outfile):
            line = line.decode()
            line = line.strip().split("\t")
            kmer, count = line[0], int(line[1])
            if count < self.histo_valley_freq:
                continue
            kmer_dictionary[kmer] = count

        return kmer_outfile, kmer_dictionary

    def __get_number_of_reads_to_analyze(self, expected_genome_size=5000000, coverage_cap=75):
        """
        This function estimates the number of reads to cover a genome of size (x) by coverage (y)
        using the median read length of the sample. It returns a tuple with the number of reads
        and whether the sample requires subsampling,

        :param expected_genome_size: The expected geneome size (in this instance a high estimate)
        :param coverage_cap: The amount of coverage allowed
        :return: number of reads to analyze: int
        """
        number_of_reads_for_75x_coverage = (expected_genome_size * coverage_cap) / self.bc_mean_read_length
        return number_of_reads_for_75x_coverage, True

    @timeit
    def __generate_kmer_histo(self):
        """
        system call to jellyfish histo to create histo output

        :return: None
        """
        if os.path.isfile(self.jellyfish_file_path) is False:
            raise ValueError("Could not find jellyfish file to generate histo values "
                             "[{0}]".format(self.jellyfish_file_path))
        histo_outfile = os.path.join(self.strain_directory, self.name + ".histo")
        try:
            check_call([self.jellyfish_app, "histo", "--output", histo_outfile, "--threads", "4",
                       self.jellyfish_file_path])
        except CalledProcessError as err:
            print(err)
        self.jellyfish_histo_path = histo_outfile
        self.kmer_histo = self.__parse_histogram()
        return

    def __parse_histogram(self):
        histo = []
        for line in open(self.jellyfish_histo_path):
            line = line.strip().split(" ")
            freq, count = int(line[0]), int(line[1])
            histo.append((freq, count,))
        return histo

    @timeit
    def count_kmers(self, bam_file=True, gzipped=False, threads=15, max_mem='1G'):
        """
        Counts the kmers using a system call to jellyfish

        :param bam_file: True if bamfile
        :param gzipped: True if gzipped file
        :param threads: number of cpus for jellyfish to call
        :param max_mem: memory size to allocate to jellyfish
        :return: path_to_jellyfish output file
        """
        jellyfish_outfile = os.path.join(self.strain_directory, self.name + ".jf")

        if bam_file:
            # see if we need to subsample

            reads_to_analyze = self.COVERAGE_CAP_READS
            # if self.COVERAGE_CAP_READS is not None and self.bc_read_count is not None:

            bam_to_fastq_process = Popen(["bamToFastq", "-i", self.bam_filepath, "-fq", "/dev/stdout"], stdout=PIPE)
            head_process = Popen(["head", "--lines", str(reads_to_analyze * 4)], stdin=bam_to_fastq_process.stdout,
                                 stdout=PIPE)
            jellyfish_process = Popen([self.jellyfish_app, "count", "--mer-len", str(self.KMER_SIZE), "--size",
                                       max_mem, "-t", str(threads), "--canonical", "--output", jellyfish_outfile,
                                       "/dev/fd/0"], stdin=head_process.stdout)

            bam_to_fastq_process.stdout.close()
            head_process.stdout.close()
            out, err = jellyfish_process.communicate()
            if out is not None:
                print(out)  # pragma: no cover
            if err is not None:  # looks good  # pragma: no cover
                print("JELLYFISH ERROR:", err)

        # NOT COVERING THIS CODE SINCE IT IS NOT NEEDED FOR PLUGIN
        else:  # pragma: no cover
            if gzipped:  # pragma: no cover
                gzip_process = Popen(["gzip", "--decompress", "--to-stdout",  self.bam_filepath], stdout=PIPE)
                jellyfish_process = Popen([self.jellyfish_app, "count", "--mer-len", str(self.KMER_SIZE),
                                           "--size", max_mem, "-t", str(threads), "--canonical", "--output",
                                           jellyfish_outfile, "/dev/fd/0"], stdin=gzip_process.stdout)

                gzip_process.stdout.close()
                out, err = jellyfish_process.communicate()
                if out is not None:
                    print(out)
                if err is not None:  # looks good
                    print("JELLYFISH ERROR:", err)

            else:  # pragma: no cover
                jellyfish_process = Popen([self.jellyfish_app, "count", "--mer-len", str(self.KMER_SIZE),
                                           "--size", max_mem, "-t", str(threads), "--canonical", "--output",
                                           jellyfish_outfile, self.bam_filepath],)
                out, err = jellyfish_process.communicate()
                if out is not None:
                    print(out)
                if err is not None:  # looks good
                    print("JELLYFISH ERROR:", err)

        return jellyfish_outfile

    @timeit
    def __get_initial_shallow(self):
        """
        finds the shallowest point between two distrubutions

        :return: index in histo, frequency, count
        """
        shallow = 1000000000000000000
        count_down = None
        frequency = 0
        index = 0
        if self.is_reference:
            return 0, 0, 0
        for i, values in enumerate(self.kmer_histo):
            freq, count = values
            if count < shallow:
                shallow = count
                index = i
                frequency = freq
                count_down = 10
            else:
                if count_down > 0:
                    count_down -= 1
                else:
                    return index, frequency, shallow
        return 1, 0, 0

    def add_basecalling_metadata(self, base_caller_directory):
        """
        adds metadata from the base calling results

         u'Q20': {
               u'max_read_length': 517,
               u'mean_read_length': 192,
               u'num_bases': 295975728,
               u'num_reads': 1539447,
            }
         u'full': {
               u'max_read_length': 535,
               u'mean_read_length': 250,
               u'num_bases': 385005746,
               u'num_reads': 1539447,
               }

        :param base_caller_directory:
        :return: (updates properties)
        """
        file_extension = "_rawlib.ionstats_basecaller.json"
        base_caller_metadata_path = os.path.join(str(base_caller_directory), self.name + file_extension)
        if os.path.isfile(base_caller_metadata_path) is False:
            raise ValueError("could not find basecalling metadata at [{0}]".format(base_caller_metadata_path))

        basecalling_metadata = json.load(open(base_caller_metadata_path))

        self.bc_max_read_length = int(basecalling_metadata['full']['max_read_length'])
        self.bc_mean_read_length = int(basecalling_metadata['full']['mean_read_length'])
        self.bc_base_count = int(basecalling_metadata['full']['num_bases'])
        self.bc_read_count = int(basecalling_metadata['full']['num_reads'])

        self.bc_max_read_length_q20 = int(basecalling_metadata['Q20']['max_read_length'])
        self.bc_mean_read_length_q20 = int(basecalling_metadata['Q20']['mean_read_length'])
        self.bc_read_count_q20 = int(basecalling_metadata['Q20']['num_reads'])

        if self.basecaller_stats is None:
            return

        self.bc_base_count_q20 = self.basecaller_stats.barcode_stats[self.barcode_name]["Q20_bases"]
        self.COVERAGE_CAP_READS, _ = self.__get_number_of_reads_to_analyze()

    def __parse_metadata(self):
        """
        parses the metadata from the sample
        :return:
        """
        if self.metadata is None:
            return

        for k, v in self.metadata.items():
            if k == 'read_count':
                self.read_count = v

            if k == 'barcode_adapter':
                self.barcode_adapter = v

            if k == 'sample':
                self.sample = v.replace(" ", "_")
                self.strain_id = v.replace(" ", "_")

            if k == 'bam_filepath':
                self.bam_filepath = v
                self.bam_file = v

            if k == 'barcode_name':
                self.barcode_name = v

            if k == 'barcode_sequence':
                self.barcode_sequence = v

            if k == 'barcode_index':
                self.barcode_index = v

            if k == 'sample_id':
                self.sample_id = v.replace(" ", "_")
                if v.strip() == "":
                    v = "control"
                self.accession = v.replace(" ", "_")
                if v.lower() == "control":
                    self.control_sample = True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Properties
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    @property
    def sample_id(self):  # THIS IS THE ACCESSION NUMBER
        return self._sample_id

    @sample_id.setter
    def sample_id(self, value):
        # TODO SHOULD WE CHECK WHETHER THE SAMPLE_ID ID VALID
        self._sample_id = value

    @property
    def sample(self):
        return self._sample

    @sample.setter
    def sample(self, value):
        self._sample = value

    @property
    def accession(self):
        return self._accession

    @accession.setter
    def accession(self, value):
        self._accession = value.replace("-", "")

    @property
    def strain_id(self):
        return self._strain_id

    @strain_id.setter
    def strain_id(self, value):
        self._strain_id = value.replace(" ", "_")

    @property
    def barcode_index(self):
        return self._barcode_index

    @barcode_index.setter
    def barcode_index(self, value):
        self._barcode_index = value

    @property
    def barcode_sequence(self):
        return self._barcode_sequence

    @barcode_sequence.setter
    def barcode_sequence(self, value):
        self._barcode_sequence = value

    @property
    def barcode_name(self):
        return self._barcode_name

    @barcode_name.setter
    def barcode_name(self, value):
        self._barcode_name = value

    @property
    def bam_file(self):
        return self._bam_file

    @bam_file.setter
    def bam_file(self, value):
        value = str(value)
        self._bam_file = os.path.basename(value)

    @property
    def bam_filepath(self):
        return self._bam_filepath

    @bam_filepath.setter
    def bam_filepath(self, value):
        value = str(value)
        if os.path.isfile(value) is False:
            raise ValueError("The bam_filepath is not valid [{0}]".format(value))
        self._bam_filepath = value

    @property
    def strain_directory(self):
        return self._strain_directory

    @strain_directory.setter
    def strain_directory(self, value):
        value = str(value)
        if os.path.isdir(value) is False:
            os.mkdir(value)
        self._strain_directory = value

    @property
    def read_count(self):
        return self._read_count

    @read_count.setter
    def read_count(self, value):
        self._read_count = int(value)

    @property
    def barcode_adapter(self):
        return self._barcode_adapter

    @barcode_adapter.setter
    def barcode_adapter(self, value):
        self._barcode_adapter = value

    @property
    def kmers(self):
        return self._kmers

    @kmers.setter
    def kmers(self, value):
        if type(value) != dict:
            raise ValueError("invalid type for kmers [{0}]".format(type(value)))
        self._kmers = value

    @property
    def coverage(self):
        return self._coverage

    @coverage.setter
    def coverage(self, value):
        if type(value) != int or value < 0:
            raise ValueError("the coverage value is invalid [{0}]".format(value))
        self._coverage = value

    @property
    def distinct_kmers(self):
        return self._distinct_kmers

    @distinct_kmers.setter
    def distinct_kmers(self, value):
        if type(value) != int or value < 0:
            raise ValueError("the distinct kmers value is invalid [{0}]".format(value))
        self._distinct_kmers = value

    @property
    def kmer_count(self):
        return self._kmer_count

    @kmer_count.setter
    def kmer_count(self, value):
        if type(value) != int or value < 0:
            raise ValueError("the kmer count is invalid [{0}]".format(value))
        self._kmer_count = value

    @property
    def histo_valley_index(self):
        return self._histo_valley_index

    @histo_valley_index.setter
    def histo_valley_index(self, value):
        if type(value) != int:
            raise ValueError("the index of the histo valley is invalid [{0}]".format(value))
        self._histo_valley_index = value

    @property
    def histo_valley_freq(self):
        return self._histo_valley_freq

    @histo_valley_freq.setter
    def histo_valley_freq(self, value):
        if type(value) != int:
            raise ValueError("the index of the histo freq is invalid [{0}]".format(value))
        self._histo_valley_freq = value

    @property
    def histo_valley_count(self):
        return self._histo_valley_count

    @histo_valley_count.setter
    def histo_valley_count(self, value):
        if type(value) != int:
            raise ValueError("the index of the histo count is invalid [{0}]".format(value))
        self._histo_valley_count = value

    @property
    def jellyfish_histo_path(self):
        return self._jellyfish_histo_path

    @jellyfish_histo_path.setter
    def jellyfish_histo_path(self, file_path):
        if os.path.isfile(file_path) is False:
            raise ValueError("jellyfish histo output file missing [{0}]".format(file_path))
        self._jellyfish_histo_path = file_path

    @property
    def kmer_path(self):
        return self._kmer_path

    @kmer_path.setter
    def kmer_path(self, file_path):
        if os.path.isfile(file_path) is False:
            raise ValueError("kmer file path is missing [{0}]".format(file_path))
        self._kmer_path = file_path

    @property
    def jellyfish_file_path(self):
        return self._jellyfish_file_path

    @jellyfish_file_path.setter
    def jellyfish_file_path(self, file_path):
        if os.path.isfile(file_path) is False:
            raise ValueError("jellyfish output file missing [{0}]".format(file_path))
        self._jellyfish_file_path = file_path

    @property
    def full_length_ssu_path(self):
        return self._full_length_ssu_path

    @full_length_ssu_path.setter
    def full_length_ssu_path(self, file_path):
        if os.path.isfile(file_path) is False:
            raise ValueError("full 16S rrna output file missing [{0}]".format(file_path))
        self._full_length_ssu_path = file_path

    @property
    def rrna_reads_path(self):
        return self._rrna_reads_path

    @rrna_reads_path.setter
    def rrna_reads_path(self, file_path):
        if os.path.isfile(file_path) is False:
            raise ValueError("16S rrna read output file missing [{0}]".format(file_path))
        self._rrna_reads_path = file_path

    @property
    def rrna_contig_path(self):
        return self._rrna_contig_path

    @rrna_contig_path.setter
    def rrna_contig_path(self, file_path):
        if os.path.isfile(file_path) is False:
            print("spades contig file output file missing [{0}]".format(file_path))
            return
        self._rrna_contig_path = file_path

    @property
    def hits(self):
        return self._hits

    @hits.setter
    def hits(self, value):
        if len(value) > 0:
            self.species_predicted_genome_size = value[0].median_genome_size_species
            self.genus_predicted_genome_size = value[0].median_genome_size_genus
        self._hits = value

    @property
    def species_predicted_genome_size(self):
        return self._species_predicted_genome_size

    @species_predicted_genome_size.setter
    def species_predicted_genome_size(self, value):
        if value is None:
            return
        if value < 0:
            raise ValueError("Invalid species genome size [{0}]".format(value))

        self._species_predicted_genome_size = value

    @property
    def genus_predicted_genome_size(self):
        return self._genus_predicted_genome_size

    @genus_predicted_genome_size.setter
    def genus_predicted_genome_size(self, value):
        if value is None:
            return
        if value < 0:
            raise ValueError("Invalid genus genome size [{0}]".format(value))
        self._genus_predicted_genome_size = value

    @property
    def kmer_backup_path(self):
        if self.backup_directory is None:
            return None
        else:
            return os.path.join(self.backup_directory, os.path.basename(self.kmer_path))

    @property
    def strain_path(self):
        return os.path.join(self.strain_directory, self.name + "_strain.json")

    @property
    def strain_backup_path(self):
        if self.backup_directory is None:
            return None
        else:
            return os.path.join(self.backup_directory, os.path.basename(self.strain_json_path))
