import os
import json
from unittest import TestCase
from strain_typing import Strain
from strain_typing import SpeciesClassifier
from strain_typing.utility_functions import get_kmers
import sys


def load_json(json_file):
    with open(json_file, 'r') as fh:
        return json.load(fh)


class TestClassifier(TestCase):

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
        strain_name = "IonCode_0101"
        metadata = load_json(os.path.join(cls.resource_dir, "barcodes.json"))[strain_name]
        startplugin = load_json(os.path.join(cls.resource_dir, "startplugin.json"))
        plugin_dir = os.path.join("/results/plugins/",  plugin_enviroment)
        cls.test_strain = Strain(strain_name, 1, False, plugin_dir,
                                 metadata=metadata, startplugin_metadata=startplugin)

        cls.classifier = SpeciesClassifier(plugin_dir, load_test_subset=True)
        cls.test_strain.strain_directory = None

    def test_load_references(self):
        self.assertEquals(5, len(self.classifier.references))

    def test_load_genome_size_genus(self):
        self.assertIn(2786690, self.classifier.genus_genome_sizes['Staphylococcus'])
        self.assertEquals(3, len(self.classifier.genus_genome_sizes.keys()))

    def test_load_genome_size_species(self):
        self.assertIn(2491640, self.classifier.species_genome_sizes['Staphylococcus lugdunensis'])
        self.assertIn('Staphylococcus lugdunensis', self.classifier.species_genome_sizes.keys())


    def test_compare_reference_classifier_hit(self):
        ref = """AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGT
                 AGCACAGAAGAGCTTGCTCTTCGGGTGACGAGCGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGAT
                 GGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCTTCGGACCAAAGAGGGGGACCTTC
                 GGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGAGGTAATGGCTCACCTAGGCGAC
                 GATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAG
                 GCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCT
                 TCGGGTTGTAAAGCACTTTCAGCGGGGAGGAAGGCGATAAGGTTAATAACCTTGTCGATTGACGTTACCC
                 GCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAA
                 TGACTGGGCGTAAAGCGCACGCAGGCGGTCTGTTAAGTTGGATGTGAAATCCCCGGGCTTAACCTGGGAA
                 CTGCATTCAAAACTGACAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCG
                 TAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGC
                 GTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGTCGACTTGGAGGTTGTGCC
                 CTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAAC
                 TCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACC
                 TTACCTACTCTTGACATCCACAGAATTCGGTAGAGATACCTTAGTGCCTTCGGGAACTGTGAGACAGGTG
                 CTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCC
                 TTTGTTGCCAGCGAGTAATGTCGGGAACTCAAAGGAGACTGCCGGTGATAAACCGGAGGAAGGTGGGGAT
                 GACGTCAAGTCATCATGGCCCTTACGAGTAGGGCTACACACGTGCTACAATGGCGTATACAAAGAGAAGC
                 GAACTCGCGAGAGCAAGCGGACCTCATAAAGTACGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCA
                 TGAAGTCGGAATCGCTAGTAATCGTAGATCAGAATGCTACGGTGAATACGTTCCCGGGCCTTGTACACAC
                 CGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCAC
                 TTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCC
                 TTA""".replace(" ", "").replace("\n", "")

        kmer_dictionary = {kmer: 10 for kmer in get_kmers(ref)}

        hit = self.classifier.compare_references(kmer_dictionary)[0]
        self.assertEquals(100, hit.kmer_identity)
        self.assertEquals("Pectobacterium parmentieri", hit.species)
        self.assertEquals("Pectobacterium", hit.genus)
        self.assertEquals("NR_153752.1", hit.accession)
        self.assertIsNone(hit.median_genome_size_species)
        self.assertIsNone(hit.median_genome_size_genus)
        self.assertEquals({10}, set(hit.counts))

    def test_compare_reference_classifier__no_hits(self):
        ref = """AATTTCCCGGGAAATCTCTCTCTCTCTCTCT""".replace(" ", "").replace("\n", "")

        kmer_dictionary = {kmer: 10 for kmer in get_kmers(ref)}

        self.assertEquals([], self.classifier.compare_references(kmer_dictionary, ))


