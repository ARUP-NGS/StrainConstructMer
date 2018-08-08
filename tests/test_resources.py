import os
import unittest
from unittest import TestCase
from strain_typing import Strain
from subprocess import Popen
from subprocess import PIPE
from hashlib import md5
from strain_typing.utility_functions import check_mongo_collection_for_record
import tempfile


def md5_checksum(file_path):
    with open(file_path, 'rb') as fh:
        m = md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()


class TestResources(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.s = Strain.__new__(Strain)  # create a class without calling __init__(); want to test class variables
        plugin_enviroment = None
        for _dir in __file__.split("/"):
            if "StrainConstructMer" in _dir:
                plugin_enviroment = _dir
                break
        if plugin_enviroment is None:
            raise OSError("could not resolve plugin and enviroment from [{0}]".format(__file__))

        cls.resource_dir = "/results/plugins/" + plugin_enviroment + "/tests/test_resources"
        cls.plugin_dir = "/results/plugins/" + plugin_enviroment

    def test_jellyfish_version(self):
        """
        test that the expected strain_instance of jellyfish is correct
        """
        expected_result = 'jellyfish 2.2.7'
        actual_results = None

        version_out = Popen([self.s.jellyfish_app, "--version"], stderr=PIPE, stdout=PIPE)  # capture output

        for line in version_out.stdout:
            if line.strip() == "":
                continue
            actual_results = line.strip()
        self.assertEqual(expected_result, actual_results)


    # def test_conda_version(self):
    #     expected_result = 'conda 4.5'
    #     actual_results = None
    #
    #     version_out = Popen(['/results/arup_pkgs/conda/anaconda3/bin/conda', "--version"], stderr=PIPE, stdout=PIPE)
    #
    #     for line in version_out.stdout:
    #         if line.strip() == "":
    #             continue
    #         actual_results = line.strip()
    #     self.assertIn(expected_result, actual_results)

    def test_barrnap_version(self):
        expected_result = 'barrnap 0.9'
        actual_results = None

        version_out = Popen([self.s.barrnap_app, "--version"], stderr=PIPE, stdout=PIPE)

        for line in version_out.stderr:
            if line.strip() == "":
                continue
            actual_results = line.strip()
        self.assertEqual(expected_result, actual_results)

    def test_sequester_version(self):
        expected_result = 'sequester2 0.9.1'
        actual_results = None

        version_out = Popen([self.s.sequester_app, "version"], stderr=PIPE, stdout=PIPE)

        for line in version_out.stderr:
            if line.strip() == "":
                continue
            actual_results = line.strip()
        self.assertEqual(expected_result, actual_results)

    def test_sequester_ssu_db(self):
        expected_result = 'd2f4742239d731630ff18545e6162f44'
        actual_results = md5_checksum(self.s.sequester_16S_db)
        self.assertEqual(expected_result, actual_results)

    def test_anaconca_python_version(self):
        expected_result = 'Python 3.6'
        version_out = Popen(['/results/arup_pkgs/conda/anaconda3/bin/python', "--version"], stderr=PIPE, stdout=PIPE)

        for line in version_out.stderr:
            if line.strip() == "":
                continue
            actual_results = line.strip()
            self.assertIn(expected_result, actual_results)

    def test_pymongo_version(self):
        expected_result = 'pymongo 3.4.0 py36_0'
        actual_results = None

        version_out = Popen(['/results/arup_pkgs/conda/anaconda3/bin/conda', "list", 'pymongo'], stderr=PIPE,
                            stdout=PIPE)

        for line in version_out.stdout:
            if line.strip() in [""]:
                continue
            if line[0] == "#":
                continue
            actual_results = " ".join(line.split())
        self.assertEqual(expected_result, actual_results)

    def test_if_database_is_reachable_dev(self):
        expected_result = (True, -33)
        md = {'sample': 'test_connection',
              'sample_id': 'test_connection',
              'bam_filepath': 'test_connection',
              'read_count': "0"}

        actual_result = check_mongo_collection_for_record(self.plugin_dir, "Strains_DEV", md, create_placeholder=False)
        self.assertEqual(expected_result, actual_result)

    def test_if_database_is_reachable_cert(self):
        expected_result = (True, -33)
        md = {'sample': 'test_connection',
              'sample_id': 'test_connection',
              'bam_filepath': 'test_connection',
              'read_count': "0"}

        actual_result = check_mongo_collection_for_record(self.plugin_dir, "Strains_CERT", md, create_placeholder=False)
        self.assertEqual(expected_result, actual_result)

    def test_if_database_is_reachable_prod(self):
        expected_result = (True, -33)
        md = {'sample': 'test_connection',
              'sample_id': 'test_connection',
              'bam_filepath': 'test_connection',
              'read_count': "0"}

        actual_result = check_mongo_collection_for_record(self.plugin_dir, "Strains_PROD", md, create_placeholder=False)
        self.assertEqual(expected_result, actual_result)

    def test_if_shared_drive_is_available_dev(self):
        try:
            tf = tempfile.NamedTemporaryFile(dir="/s5xl_shared/Strains_DEV/", delete=True)
            actual_result = os.path.exists(tf.name)
        except OSError:
            actual_result = False

        self.assertTrue(actual_result)

    def test_if_shared_drive_is_available_cert(self):
        try:
            tf = tempfile.NamedTemporaryFile(dir="/s5xl_shared/Strains_CERT/", delete=True)
            actual_result = os.path.exists(tf.name)
        except OSError:
            actual_result = False

        self.assertTrue(actual_result)

    def test_if_shared_drive_is_available_prod(self):
        try:
            tf = tempfile.NamedTemporaryFile(dir="/s5xl_shared/Strains_PROD/", delete=True)
            actual_result = os.path.exists(tf.name)
        except OSError:
            actual_result = False

        self.assertTrue(actual_result)
