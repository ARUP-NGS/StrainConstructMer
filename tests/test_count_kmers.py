from unittest import TestCase
from strain_typing import Strain
import tempfile
import shutil
import os
from hashlib import md5


def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()


class TestCountKmers(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_strain = Strain.__new__(Strain)
        cls.test_strain.strain_directory = tempfile.mkdtemp()
        cls.test_strain.name = "test_count"

        plugin_enviroment = None
        for _dir in __file__.split("/"):
            if "StrainConstructMer" in _dir:
                plugin_enviroment = _dir
                break
        if plugin_enviroment is None:
            raise OSError("could not resolve plugin and enviroment from [{0}]".format(__file__))

        cls.plugin_name = plugin_enviroment
        cls.resource_dir = "/results/plugins/" + plugin_enviroment + "/tests/test_resources"

    def test_counting_kmers_with_bamfile(self):
        self.test_strain.bam_filepath = os.path.join("/results/plugins/", self.plugin_name,
                                                     "tests/test_resources/lambda_NC_001416_simulation.bam")
        jf_file = self.test_strain.count_kmers(bam_file=True, threads=5, max_mem="50M")
        self.assertTrue(os.path.isfile(jf_file))

    def test_counting_kmers_with_malformed_bamfile(self):
        self.test_strain.bam_filepath = os.path.join("/results/plugins/", self.plugin_name,
                                                     "tests/test_resources/lambda_NC_001416_simulation_malformed.bam")
        jf_file = self.test_strain.count_kmers(bam_file=True, threads=5, max_mem="50M")
        self.assertTrue(os.path.isfile(jf_file))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.test_strain.strain_directory)
