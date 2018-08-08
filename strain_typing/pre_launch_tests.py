import os
import tempfile
from subprocess import Popen
from subprocess import PIPE
from .utility_functions import check_mongo_collection_for_record
from .strain import Strain


def test_jellyfish_version():
    """
    Sanity check to make sure that jellyfish is installed and that it is the correct version.

    :return:  passed(bool), expected results(str), actual result(str)
    """
    s = Strain.__new__(Strain)
    expected_result = 'jellyfish 2.2.7'
    actual_results = None

    version_out = Popen([s.jellyfish_app, "--version"], stderr=PIPE, stdout=PIPE)  # capture output

    for line in version_out.stdout:
        if line.strip() == "":
            continue
        actual_results = line.strip()
    return actual_results == expected_result, expected_result, actual_results, "jellyfish_version"


def test_anaconca_python_version():
    """
    Sanity check to make sure that anaconda python is right version.

    :return:  passed(bool), expected results(str), actual result(str)
    """
    expected_result = 'Python 3.6'
    actual_results = None

    version_out = Popen(['/results/arup_pkgs/conda/anaconda3/bin/python', "--version"], stderr=PIPE, stdout=PIPE)

    for line in version_out.stderr:
        if line.strip() == "":
            continue
        actual_results = line.strip()
    return expected_result in actual_results, expected_result, actual_results, "python_version"


def test_pymongo_version():
    expected_result = 'pymongo 3.4.0 py36_0'
    actual_results = None

    version_out = Popen(['/results/arup_pkgs/conda/anaconda3/bin/conda', "list", 'pymongo'], stderr=PIPE, stdout=PIPE)

    for line in version_out.stdout:
        if line.strip() in [""]:
            continue
        if line[0] == "#":
            continue
        actual_results = " ".join(line.split())
    return actual_results == expected_result, expected_result, actual_results, "pymongo version"


def test_if_database_is_reachable(plugin_directory, collection):
    expected_result = (True, -33)
    md = {'sample': 'test_connection',
          'sample_id': 'test_connection',
          'bam_filepath': 'test_connection',
          'read_count': "0"}

    actual_result = check_mongo_collection_for_record(plugin_directory, collection, md, create_placeholder=False)
    return actual_result == expected_result, expected_result, actual_result, "mongo_connection"


def test_if_shared_drive_is_available(backup_location):
    expected_result = True
    try:
        tf = tempfile.NamedTemporaryFile(dir=backup_location, delete=True)
        actual_result = os.path.exists(tf.name)
    except OSError:
        actual_result = False

    return actual_result == expected_result, expected_result, actual_result, "backup_network_avaliable"
