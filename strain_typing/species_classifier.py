import os
from Bio import SeqIO
from strain_typing.utility_functions import get_kmers
from strain_typing.utility_functions import timeit
import csv
from collections import defaultdict
import numpy
import sys


class ClassifierHit(object):
    def __init__(self, kmer_identity, species, genus, accession, median_genome_size_species, median_genome_size_genus,
                 min_genome_size_species, min_genome_size_genus, max_genome_size_species, max_genome_size_genus,
                 counts):

        self.kmer_identity = kmer_identity
        self.species = species
        self.genus = genus
        self.accession = accession
        self.median_genome_size_species = median_genome_size_species
        self.median_genome_size_genus = median_genome_size_genus
        self.min_genome_size_species = min_genome_size_species
        self.min_genome_size_genus = min_genome_size_genus
        self.max_genome_size_species = max_genome_size_species
        self.max_genome_size_genus = max_genome_size_genus
        self.counts = counts

    def formatted_output(self):
        _d = {}
        pass
        _d['kmer_identity'] = self.kmer_identity
        _d['species'] = self.species
        _d['genus'] = self.genus
        _d['accession'] = self.accession
        _d['median_genome_size_species'] = self.median_genome_size_species
        _d['median_genome_size_genus'] = self.median_genome_size_genus
        _d['min_genome_size_species'] = self.min_genome_size_species
        _d['min_genome_size_genus'] = self.min_genome_size_genus
        _d['max_genome_size_species'] = self.max_genome_size_species
        _d['max_genome_size_genus'] = self.max_genome_size_genus
        return _d

    def print_hit(self, delimiter=","):
        header = ["pct_kmer_identity",
                  "species",
                  "genus",
                  "accession",
                  "median_species_size",
                  "minimum_species_size",
                  "maximum_species_size",
                  "median_genus_size",
                  "minimum_genus_size",
                  "maximum_genus_size"]

        hit = ["{0:.1f}".format(self.kmer_identity),
               "{0}".format(self.species),
               "{0}".format(self.genus),
               "{0}".format(self.accession),
               "{0}".format(self.median_genome_size_species),
               "{0}".format(self.min_genome_size_species),
               "{0}".format(self.max_genome_size_species),
               "{0}".format(self.median_genome_size_genus),
               "{0}".format(self.min_genome_size_genus),
               "{0}".format(self.max_genome_size_genus)]
        return delimiter.join(header).upper(), delimiter.join(hit)

    def __repr__(self):
        s = "kID: {0:.1f}%\t".format(self.kmer_identity)
        s += "SP: {0}\t".format(self.species)
        s += "GN: {0}\t".format(self.genus)
        s += "ACC: {0}\t".format(self.accession)
        if self.median_genome_size_species is not None:
            s += "MED_SZ_SPECIES: {0:,}\t".format(self.median_genome_size_species)
        else:
            s += "MED_SZ_SPECIES: NA\t"
        if self.median_genome_size_genus is not None:
            s += "MED_SZ_GENUS: {0:,}\t".format(self.median_genome_size_genus)
        else:
            s += "MED_SZ_GENUS: NA"
        return s


class SpeciesClassifier(object):
    rel_reference_path = "strain_typing/resources/references_ssu/" \
                         "16S_bacteria_PRJNA33175.fa"
    rel_processed_path = "strain_typing/resources/references_ssu/" \
                         "16S_bacteria_PRJNA33175_processed.txt"
    rel_genome_path = "strain_typing/resources/references_ssu/prokaryotes.csv"

    def __init__(self, plugin_path, load_test_subset=False):
        self.reference_path = self.__set_reference_path(plugin_path)
        self.processed_path = self.__set_processed_path(plugin_path)
        self.genome_path = self.__set_genome_path(plugin_path)
        self.references = self.__load_references(load_test_subset)
        self.genus_genome_sizes, self.species_genome_sizes = self.__load_genome_sizes()

    @timeit
    def __load_genome_sizes(self):
        genome_sizes_species = defaultdict(list)
        genome_sizes_genus = defaultdict(list)

        for line in csv.reader(open(self.genome_path)):
            if line[0][0] == "#":
                continue
            if "sp." in line[0]:
                continue
            species = " ".join(line[0].split(" ")[0:2])
            genus = species.split(" ")[0]
            size = int(float(line[7]) * 1000000)
            genome_sizes_species[species].append(size)
            genome_sizes_genus[genus].append(size)
        return genome_sizes_genus, genome_sizes_species

    def __set_genome_path(self, plugin_path):
        reference_path = os.path.join(plugin_path, self.rel_genome_path)
        if os.path.isfile(reference_path) is False:
            raise ValueError("the is_reference file could not be found [{0}]".format(reference_path))
        return reference_path

    def __set_reference_path(self, plugin_path):
        reference_path = os.path.join(plugin_path, self.rel_reference_path)
        if os.path.isfile(reference_path) is False:
            raise ValueError("the is_reference file could not be found [{0}]".format(reference_path))
        return reference_path

    def __set_processed_path(self, plugin_path):
        reference_path = os.path.join(plugin_path, self.rel_processed_path)
        return reference_path

    @timeit
    def __load_references(self, load_test_subset=False):
        references = []
        if os.path.isfile(self.processed_path):
            sys.stderr.write("INFO: Found processed references, loading.\n")
            with open(self.processed_path) as rf:
                for line in rf:
                    accession, genus, species, sequence, kmer_list = line.split("\t")
                    kmer_list = kmer_list.split(";")
                    r = ClassifierReference.__new__(ClassifierReference)
                    r.accession = accession
                    r.genus = genus
                    r.species = species
                    r.sequence = sequence
                    r.kmer_list = kmer_list
                    references.append(r)
            return references
        sys.stderr.write("INFO: No processed reference file loading references and creating processed reference set.\n")
        for reference in SeqIO.parse(self.reference_path, "fasta"):
            r = ClassifierReference(reference)
            references.append(r)
            if load_test_subset and len(references) > 4:  # this is just in place to load a subset for testing
                break
        with open(self.processed_path, "w") as wf:
            for s in references:
                wf.write("{0}\n".format(str(s)))
        return references

    @timeit
    def compare_references(self, kmer_dictionary):
        hits = []
        for ref in self.references:
            results = ref.report_similarity(kmer_dictionary)
            if results is not None:
                kmer_identity, species, genus, accession, counts = results
                h = ClassifierHit(kmer_identity, species, genus, accession,
                                  self.get_median_species_genome_size(species),
                                  self.get_median_genus_genome_size(genus),
                                  self.get_min_species_genome_size(species),
                                  self.get_min_genus_genome_size(genus),
                                  self.get_max_species_genome_size(species),
                                  self.get_max_genus_genome_size(genus),
                                  counts)
                hits.append(h)
        return sorted(hits, key=lambda x: x.kmer_identity, reverse=True)

    def get_median_genus_genome_size(self, genus):
        if genus in self.genus_genome_sizes:
            return numpy.median(self.genus_genome_sizes[genus])
        else:
            return None

    def get_median_species_genome_size(self, species):
        if species in self.species_genome_sizes:
            return numpy.median(self.species_genome_sizes[species])
        else:
            return None

    def get_min_genus_genome_size(self, genus):
        if genus in self.genus_genome_sizes:
            return numpy.min(self.genus_genome_sizes[genus])
        else:
            return None

    def get_min_species_genome_size(self, species):
        if species in self.species_genome_sizes:
            return numpy.min(self.species_genome_sizes[species])
        else:
            return None

    def get_max_genus_genome_size(self, genus):
        if genus in self.genus_genome_sizes:
            return numpy.max(self.genus_genome_sizes[genus])
        else:
            return None

    def get_max_species_genome_size(self, species):
        if species in self.species_genome_sizes:
            return numpy.max(self.species_genome_sizes[species])
        else:
            return None


class ClassifierReference(object):

    def __init__(self, seqio_object):
        self.accession, self.genus, self.species = self.__parse_description(seqio_object)
        self.sequence = seqio_object.seq
        self.kmer_list = get_kmers(str(seqio_object.seq))

    def report_similarity(self, kmer_dictionary, cutoff=70):
        matches = 0
        counts = []
        for i, kmer in enumerate(self.kmer_list):
            if i > 30 and matches < 10:  # abandon if first 30 hits don't match > 33%
                return None
            if kmer in kmer_dictionary:
                counts.append(kmer_dictionary[kmer])
                matches += 1
            else:
                counts.append(0)

        kmer_identity = float(matches)/float(len(self)) * 100.0
        if kmer_identity < cutoff:
            return None

        return kmer_identity, self.species, self.genus, self.accession, counts

    @staticmethod
    def __parse_description(seqio_object):
        description_line = seqio_object.description.split(" ")
        accession = description_line[0]
        genus = description_line[1]
        species = "{0} {1}".format(description_line[1], description_line[2])
        return accession, genus, species

    def __repr__(self):
        s = "{0}\t{1}\t{2}\t".format(self.accession, self.genus, self.species)
        s += "{0}\t".format(self.sequence)
        s += "{0}".format(";".join(self.kmer_list))
        return s

    def __len__(self):
        return len(self.kmer_list)


if __name__ == '__main__':
    sc = SpeciesClassifier("/Users/331-SimmonkLPTP/git_repos_thermo/StrainConstructMer/")
    with open("/Users/331-SimmonkLPTP/git_repos_thermo/StrainConstructMer/strain_typing/"
              "resources/references_ssu/16S_bacteria_PRJNA33175_processed.txt", "w") as wf:
        for s in sc.references:
            wf.write("{0}\n".format(str(s)))
