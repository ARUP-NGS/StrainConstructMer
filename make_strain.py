import argparse
from strain_typing import Strain
from strain_typing import SpeciesClassifier
from strain_typing import AMRDatabase
from strain_typing import MLST
import os
import inspect
import sys
import datetime
import simplejson

# set up for Django rendering with TSS installed apps
from django.conf import settings
from django.conf import global_settings
from django.template.loader import render_to_string
global_settings.LOGGING_CONFIG = None

"""
This is the commandline interface to run create strains outside of the plugin
"""


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="a bam/fastq/fasta file, the file can be compressed")
    parser.add_argument("-gz", "--gzipped", action="store_true", default=False, help="The file is gzipped compressed")
    parser.add_argument("-t", "--threads", type=int, default=2, help="the number of threads to count kmers")
    parser.add_argument("-s", "--memory_size", default="500M", help="the size of the memory hash")
    parser.add_argument("-n", "--sample_name", help="name of the strain", required=True)
    parser.add_argument("-id", "--sample_id", help="sample id", required=True)
    parser.add_argument("-bf", "--bam_file", help="The file is a bamfile", default=False, action='store_true')
    parser.add_argument("-f", "--force", help="will overwrite existing output", default=False, action="store_true")
    return parser.parse_args()


def main():
    args = arguments()
    plugin_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    s = Strain("{0}_{1}".format(args.sample_name, args.sample_id), 1, False, plugin_path)

    s.sample = "{0}_{1}".format(args.sample_name, args.sample_id)
    s.strain_id = s.sample
    s.sample_id = args.sample_id
    s.accession = s.strain_id

    # create directory for output
    path_to_create = "StrainTypeMerStrain_{0}".format(s.sample)
    if os.path.isdir(path_to_create):
        if args.force:
            sys.stderr.write("WARNING OVERWRITING EXISTING DATA\n")
        else:
            # exit the directory already exists
            raise ValueError("the specified sample_name [{0}] and sample_id [{1}] have already been run".format(
                args.sample_name, args.sample_id
            ))
    else:
        os.mkdir(path_to_create)

    s.strain_directory = os.path.abspath(path_to_create)
    s.bam_file = os.path.basename(args.file)

    # initialize the AMRDatabase
    antibiotic_gene_database = AMRDatabase(plugin_path)

    # initialize the  SpeciesClassifier
    species_classifier = SpeciesClassifier(plugin_path)

    # initialize the MLST profiler
    mlst_profiler = MLST(plugin_path)

    if not settings.configured:
        settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
                           INSTALLED_APPS=('django.contrib.humanize',),
                           TEMPLATE_DIRS=(os.path.join(plugin_path, 'templates'),))

    s.bam_filepath = os.path.abspath(args.file)
    sys.stderr.write("counting kmers with {0} threads\n".format(args.threads))
    s.jellyfish_file_path = s.count_kmers(threads=args.threads, bam_file=args.bam_file, gzipped=args.gzipped,
                                          max_mem=args.memory_size)

    sys.stderr.write("done counting kmers\n")

    # s.generate_kmer_histo()
    s.setup_strain()
    s.antibiotic_genes = antibiotic_gene_database.find_antibiotic_genes(s)
    s.write_antibiotic_genes()
    s.mlst_profiles, profiles_loaded = mlst_profiler.find_mlst_profiles(s)
    s.hits = species_classifier.compare_references(s.kmers)
    #s.sequester_rrna_reads(bam_file=args.bam_file, gzipped=args.gzipped)
    #s.assemble_rrna_reads()
    #s.extract_rrna_sequence()
    s.write_hits()
    s.write_mlst()
    s.write_hits_to_html()
    s.write_results()
    s.create_resistant_table()
    # s.create_symlinks()



    render_context = {
        "autorefresh": False,
        "run_name": "test_table",
        "plugin_name": "StrainConstructMer_DEV",
        "run_user_name": "CLI",
        "plugin_user_name": "CLI",
        "strain_page_link" : '',
        "date_run": str(datetime.date.today()),
        "plugin_results": simplejson.dumps(plugin_path),
    }




    with open(os.path.join(s.strain_directory, "{0}_strain_summary.html".format(s.sample)), "w") as html_fp:
        html_fp.write(render_to_string("barcode_summary.html", render_context))


if __name__ == '__main__':
    main()
