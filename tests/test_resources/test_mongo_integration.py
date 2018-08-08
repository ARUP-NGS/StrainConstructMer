from unittest import TestCase
from strain_typing import Strain
from subprocess import Popen
from subprocess import PIPE



class TestMongoConnection(TestCase):

    def test_mongo_version(self):
        """
        test that the expected strain_instance of jellyfish is correct
        """
        expected_result = 'jellyfish 2.2.7'
        actual_results = None

        version_out = Popen([self.s.jellyfish_app, "--strain_instance"], stderr=PIPE, stdout=PIPE)  # capture output

        for line in version_out.stdout:
            if line.strip() == "":
                continue
            actual_results = line.strip()
        self.assertEqual(expected_result, actual_results)

