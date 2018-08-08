from unittest import TestCase
from strain_typing import Strain
import tempfile
import json
import os
import shutil


def load_json(json_file):
    with open(json_file, 'r') as fh:
        return json.load(fh)


class TestStrain(TestCase):

    # def __init__(self, name, strain_instance, already_processed, plugin_directory, metadata=None,
    #                  startplugin_metadata=None, run_info=None, is_reference=False):

    @classmethod
    def setUp(cls):
        plugin_enviroment = None
        for _dir in __file__.split("/"):
            if "StrainConstructMer" in _dir:
                plugin_enviroment = _dir
                break
        if plugin_enviroment is None:
            raise OSError("could not resolve plugin and enviroment from [{0}]".format(__file__))

        cls.resource_dir = "/results/plugins/" + plugin_enviroment + "/tests/test_resources"

        #cls.test_strain = Strain.__new__(Strain)
        strain_name = "IonCode_0101"
        metadata = load_json(os.path.join(cls.resource_dir, "barcodes.json"))[strain_name]
        startplugin = load_json(os.path.join(cls.resource_dir, "startplugin.json"))
        cls.test_strain = Strain(strain_name, 1, False, os.path.join("/results/plugins/", plugin_enviroment),
                                 metadata=metadata, startplugin_metadata=startplugin)
        cls.test_strain.strain_directory = None

    def test_class_varible_kmer_size(self):
        expected_result = 31
        self.assertEqual(expected_result, self.test_strain.KMER_SIZE)

    def test_class_varible_coverage_cutoff(self):
        expected_result = 25
        self.assertEqual(expected_result, self.test_strain.COVERAGE_CUTOFF)

    def test_class_varible_genome_size_cutoff(self):
        expected_result = 1000000  # 1 million
        self.assertEqual(expected_result, self.test_strain.GENOME_SIZE_CUTOFF)

    def test_class_varible_difference_from_expected_cutoff(self):
        expected_result = 15  # 1 million
        self.assertEqual(expected_result, self.test_strain.DIFFERENCE_FROM_EXPECTED_CUTOFF)

    def test_class_varible_pct_q20_bases_cutoff(self):
        expected_result = 70  # 1 million
        self.assertEqual(expected_result, self.test_strain.PCT_QC20_BASES_CUTOFF)


    def test_add_basecalling_metadata(self):

        self.test_strain.name = "IonCode_0101"
        self.test_strain.add_basecalling_metadata(self.resource_dir)

        # full values
        self.assertEqual(615, self.test_strain.bc_max_read_length)
        self.assertEqual(281, self.test_strain.bc_mean_read_length)
        self.assertEqual(223985064, self.test_strain.bc_base_count)
        self.assertEqual(795328, self.test_strain.bc_read_count)

        # Q20 values are not equal with full values
        self.assertNotEqual(615, self.test_strain.bc_max_read_length_q20)
        self.assertNotEqual(281, self.test_strain.bc_mean_read_length_q20)
        self.assertNotEqual(223985064, self.test_strain.bc_base_count_q20)

        # Q20
        self.assertEqual(558, self.test_strain.bc_max_read_length_q20)
        self.assertEqual(234, self.test_strain.bc_mean_read_length_q20)
        # self.assertEqual(186876513, self.test_strain.bc_base_count_q20)
        self.assertEqual(795328, self.test_strain.bc_read_count_q20)

        # Q17
        # self.assertEqual(615, self.test_strain.bc_max_read_length_q17)
        # self.assertEqual(274, self.test_strain.bc_mean_read_length_q17)
        # self.assertEqual(218235043, self.test_strain.bc_base_count_q17)
        # self.assertEqual(795328, self.test_strain.bc_read_count_q17)


    def test_metadata_parsed_correctly_read_count(self):
        expected_value = 795328
        self.assertEqual(expected_value, self.test_strain.read_count)

    def test_metadata_parsed_correctly_sample(self):
        expected_value = 'ST-A204-094'
        self.assertEqual(expected_value, self.test_strain.sample)

    def test_metadata_parsed_correctly_barcode_adapter(self):
        expected_value = "GGTGAT"
        self.assertEqual(expected_value, self.test_strain.barcode_adapter)

    def test_metadata_parsed_correctly_bam_filepath(self):
        expected_value = "/s5xl_shared/test_resource/" \
                         "IonCode_0101_rawlib.basecaller.bam"
        self.assertEqual(expected_value, self.test_strain.bam_filepath)

    def test_metadata_parsed_correctly_barcode_name(self):
        expected_value = "IonCode_0101"
        self.assertEqual(expected_value, self.test_strain.barcode_name)

    def test_metadata_parsed_correctly_barcode_sequence(self):
        expected_value = "CTAAGGTAAC"
        self.assertEqual(expected_value, self.test_strain.barcode_sequence)

    def test_metadata_parsed_correctly_barcode_index(self):
        expected_value = 1
        self.assertEqual(expected_value, self.test_strain.barcode_index)

    def test_metadata_parsed_correctly_sample_id(self):
        expected_value = "16-112-400676"
        self.assertEqual(expected_value, self.test_strain.sample_id)

    def test_add_strain_directory(self):
        self.test_strain.strain_directory = self.resource_dir
        self.assertTrue(os.path.isdir(self.test_strain.strain_directory))

    def test_add_strain_directory_not_valid(self):
        self.test_strain.strain_directory = tempfile.mkdtemp()
        self.assertTrue(os.path.isdir(self.test_strain.strain_directory))
        shutil.rmtree(self.test_strain.strain_directory)
        self.assertFalse(os.path.isdir(self.test_strain.strain_directory))

    # validate the create sym link function
    def test_create_symlink_with_empty_sample_name(self):
        self.test_strain.strain_directory = tempfile.mkdtemp()
        self.test_strain.sample = ""
        self.assertFalse(self.test_strain.create_symlinks())
        self.test_strain.sample = None
        self.assertFalse(self.test_strain.create_symlinks())

        shutil.rmtree(self.test_strain.strain_directory)
        if self.test_strain.named_strain_directory is not None:
            os.unlink(self.test_strain.named_strain_directory)

    def test_creates_symlink_directory_created(self):
        self.test_strain.strain_directory = tempfile.mkdtemp()
        self.assertIsNone(self.test_strain.named_strain_directory)
        self.test_strain.create_symlinks()
        self.assertIsNotNone(self.test_strain.named_strain_directory)
        self.assertTrue(os.path.isdir(self.test_strain.strain_directory))
        self.assertTrue(os.path.islink(self.test_strain.named_strain_directory))
        shutil.rmtree(self.test_strain.strain_directory)
        if self.test_strain.named_strain_directory is not None:
            os.unlink(self.test_strain.named_strain_directory)

    def test_create_symlink_for_all_strain_files(self):
        self.test_strain.strain_directory = tempfile.mkdtemp()
        tf = tempfile.NamedTemporaryFile(dir=self.test_strain.strain_directory, delete=False)
        tf_bam = tempfile.NamedTemporaryFile(dir=self.test_strain.strain_directory, delete=False)
        self.test_strain.rrna_reads_path = tf.name
        self.test_strain.bam_filepath = tf_bam.name
        self.test_strain.create_symlinks()
        self.assertTrue(os.path.isfile(self.test_strain.rrna_reads_path))
        prefix = "{0}__{1}__{2}".format(self.test_strain.accession, self.test_strain.strain_id,
                                        self.test_strain.strain_version).replace(" ", "_")
        linked_sample_name = os.path.join(self.test_strain.strain_directory,
                                          "{0}_{1}".format(prefix,
                                                           os.path.basename(self.test_strain.rrna_reads_path)))
        self.assertTrue(os.path.isfile(linked_sample_name))
        self.assertTrue(os.path.islink(linked_sample_name))

        shutil.rmtree(self.test_strain.strain_directory)
        if self.test_strain.named_strain_directory is not None:
            os.unlink(self.test_strain.named_strain_directory)

    def test_write_dictionary(self):
        self.test_strain.strain_directory = tempfile.mkdtemp()

        self.test_strain.write_results()
        self.assertTrue(os.path.isfile(self.test_strain.strain_json_path))
        shutil.rmtree(self.test_strain.strain_directory)
        if self.test_strain.named_strain_directory is not None:
            os.unlink(self.test_strain.named_strain_directory)

    def test_results_summary(self):
        self.test_strain.strain_directory = tempfile.mkdtemp()
        self.assertIsInstance(self.test_strain.result_summary(), dict)

    def test_qc_icon_is_alert(self):
        self.test_strain.sample = "BARCODE_CONTAMINATE"
        self.test_strain.accession = "BARCODE_CONTAMINATE"
        self.assertEqual(self.test_strain.alert_icon, self.test_strain.get_qc_icon())

    def test_strain_qc_fail_coverage(self):
        COVERAGE_CUTOFF = 25
        GENOME_SIZE_CUTOFF = 1000000
        # DIFFERENCE_FROM_EXPECTED_CUTOFF = 15
        PCT_QC20_BASES_CUTOFF = 70
        COVERAGE_CAP_READS = 10000000

        self.test_strain.coverage = 24
        self.test_strain.distinct_kmers = 100000000
        self.assertFalse(self.test_strain.strain_passed_qc())

    def test_strain_qc_pass_coverage(self):
        self.test_strain.coverage = 25
        self.test_strain.distinct_kmers = 100000000
        self.assertTrue(self.test_strain.strain_passed_qc())

    def test_strain_qc_pass_genome_size(self):
        self.test_strain.coverage = 25
        self.test_strain.distinct_kmers = 100000000
        self.assertTrue(self.test_strain.strain_passed_qc())

    def test_strain_qc_fail_genome_size(self):
        self.test_strain.coverage = 25
        self.test_strain.distinct_kmers = 999999
        self.assertFalse(self.test_strain.strain_passed_qc())

    def test_strain_within_expected_range(self):
        class hit_fake:
            def __init__(self):
                self.min_genome_size_species = 4000000
                self.max_genome_size_species = 5000000
                self.median_genome_size_species = 4500000
                self.min_genome_size_genus = 4000000
                self.max_genome_size_genus = 5000000
                self.median_genome_size_genus = 4500000

        self.test_strain.hits = [hit_fake()]
        self.test_strain.coverage = 25

        self.test_strain.distinct_kmers = 4000000
        self.assertTrue(self.test_strain.strain_passed_qc())

        self.test_strain.distinct_kmers = 3600000
        self.assertTrue(self.test_strain.strain_passed_qc())

        self.test_strain.distinct_kmers = 3599999
        self.assertFalse(self.test_strain.strain_passed_qc())

        self.test_strain.distinct_kmers = 5000000
        self.assertTrue(self.test_strain.strain_passed_qc())

        self.test_strain.distinct_kmers = 5500000
        self.assertTrue(self.test_strain.strain_passed_qc())

        self.test_strain.distinct_kmers = 5500001
        self.assertFalse(self.test_strain.strain_passed_qc())
