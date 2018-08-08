import sys
import time
import numpy
from subprocess import Popen, PIPE
from itertools import chain
import os
import sqlite3


def find_peaks(data, spacing=1, limit=None):
    """Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return: peak indices
    """
    length = data.size
    x = numpy.zeros(length + 2 * spacing)
    x[: spacing] = data[0] - 1.e-6
    x[-spacing:] = data[-1] - 1.e-6
    x[spacing: spacing + length] = data
    peak_candidate = numpy.zeros(length)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start: start + length]  # before
        start = spacing
        h_c = x[start: start + length]  # central
        start = spacing + s + 1
        h_a = x[start: start + length]  # after
        peak_candidate = numpy.logical_and(peak_candidate, numpy.logical_and(h_c > h_b, h_c > h_a))

    ind = numpy.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
    if limit is not None:
        ind = ind[data[ind] > limit]
    return ind


def create_sql_database(db_path):
    """
    creates the sqllite database that holds information about the strain to perform comparison

    :param db_path: the full path (including name to the database)
    :return: True if created false otherwise
    """
    if os.path.exists(db_path):
        return False

    db = sqlite3.connect(db_path)
    cur = db.cursor()
    cur.execute("""CREATE TABLE strains(sample TEXT NOT NULL,
                                        sample_id TEXT,
                                        sample_instance INT NOT NULL,
                                        read_count INT NOT NULL,
                                        bam_file_path TEXT NOT NULL,
                                        passed_qc BOOLEAN NOT NULL,
                                        kmer_file_path TEXT NOT NULL,
                                        kmer_backup_path TEXT,
                                        strain_file_path TEXT NOT NULL,
                                        strain_backup_path TEXT,
                                        PRIMARY KEY (sample, sample_id, sample_instance),
                                        CONSTRAINT composite_name UNIQUE (sample, sample_id, read_count, bam_file_path)
                                        );""")

    cur.execute("CREATE INDEX index_read_count ON strains (read_count);")
    cur.execute("CREATE INDEX index_bam_file_path ON strains (bam_file_path);")

    db.commit()
    db.close()
    if os.path.exists(db_path):
        return True
    return False


def check_sql_db_for_record(db_path, sample, sample_id, read_count, bam_file_path):
    """
    This checks to see if a record already exists for a sample using the sample, sample_id, read_count,
    and bam_file_path.

    (sample, sample_id, read_count, bam_file_path) should be universally unique

    :param db_path: The sqllite database path
    :param sample: sample from template
    :param sample_id: sample_id from template
    :param read_count: read count of the sample
    :param bam_file_path: bam file path for the sample
    :return: record_exists (bool), strain_instance (int,None)
    """
    db = sqlite3.connect(db_path)
    cur = db.cursor()
    cur.execute("""SELECT * FROM strains WHERE sample = (?) 
                                           AND sample_id = (?) 
                                           AND read_count = (?) 
                                           AND bam_file_path = (?);""", (sample, sample_id, read_count, bam_file_path,))
    data = cur.fetchall()

    # record exist return instance with same read count and bam_file_path
    if len(data) > 0:
        # need to check that the files exist [plugin results may have been deleted]
        sample, sample_id, sample_instance, read_count, bam_file_path, passed_qc, kmer_file_path, kmer_backup_path, \
        strain_file_path, strain_backup_path = data[0]

        if os.path.isfile(kmer_file_path) is False and os.path.isfile(strain_file_path) is False:
            # the strain record has been deleted need to remove record and replace
            cur.execute("""DELETE FROM strains WHERE sample = (?) 
                                                 AND sample_id = (?) 
                                                 AND read_count = (?) 
                                                 AND bam_file_path = (?);""",
                        (sample, sample_id, read_count, bam_file_path,))
            db.commit()
            db.close()
            return False, sample_instance
        db.close()
        return True, sample_instance

    # record does not exist see if it is a new instance
    cur.execute("""SELECT sample_instance FROM strains WHERE sample = (?) AND sample_id = (?)
                   ORDER BY sample_instance DESC;""", (sample, sample_id,))

    data = cur.fetchall()
    if len(data) == 0:
        db.close()
        return False, 1
    db.close()
    return False, data[0][0] + 1

def add_record_to_sql(db_path, strain):
    """
    Inserts the following info into the sqllite database

    sample
    sample_id
    sample_instance
    read_count
    bam_file_path
    passed_qc
    kmer_file_path
    kmer_backup_path
    strain_file_path
    strain_backup_path

    :param dbpath:
    :param strain:
    :return:
    """
    try:
        db = sqlite3.connect(db_path)
        cur = db.cursor()
        cur.execute("""INSERT OR REPLACE INTO strains (sample, sample_id, sample_instance, read_count, bam_file_path, 
                                                       passed_qc, kmer_file_path, kmer_backup_path, strain_file_path, 
                                                       strain_backup_path) 
                                            VALUES (?,?,?,?,?,?,?,?,?,?);""", (strain.sample,
                                                                               strain.sample_id,
                                                                               strain.strain_instance,
                                                                               strain.read_count,
                                                                               strain.bam_filepath,
                                                                               strain.strain_passed_qc(),
                                                                               strain.kmer_path,
                                                                               strain.kmer_backup_path,
                                                                               strain.strain_json_path,
                                                                               strain.strain_backup_path))
        db.commit()
        db.close()
        return "SUCCESS"
    except:
        return "FAILED"


def timeit(method):
    """
    prints out the execution time to sys.stderr of the decorated function
    :param method: function to decorate
    :return: decorated function
    """

    def timed(*args, **kw):
        start = time.time()
        result = method(*args, **kw)
        end = time.time()
        time.sleep(0.1)
        sys.stderr.write("INFO: {0:.2f} seconds to run '{1}' function\n".format(end - start, method.__name__))
        return result
    return timed


def reverse_complement(sequence):
    """
    return the reverse complement of the sequence string
    :param sequence:
    :return:
    """

    base_map = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "Y": "R", "R": "Y", "M":
                "K", "K": "M", "S": "W", "W": "S", "V": "D", "D": "V", "H": "B", "B": "H"}

    try:
        new_sequence = "".join([base_map[char] for char in reversed(sequence)])
    except KeyError as e:
        raise KeyError("{0} char not present in ".format(e))
    return new_sequence


def get_kmers(sequence, kmer_size=31, canonical_kmers=True):
    """
    Creates a list of kmer found in the sequence

    :param sequence: String of the sequence
    :param kmer_size: the length of the kmer
    :param canonical_kmers: Will return the kmer orientation that is the smallest
    :return: set of kmers in the sequence
    """
    sequence = sequence.upper()
    kmers = []
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]
        if kmer.count("A") + kmer.count("T") + kmer.count("G") + kmer.count("C") != len(kmer):
            # IGNORE KMER
            continue
        if canonical_kmers:
            rckmer = reverse_complement(kmer)
            if rckmer < kmer:
                kmers.append(rckmer)
            else:
                kmers.append(kmer)
        else:
            kmers.append(kmer)
    return kmers
