import os
try:
    import cPickle as pickle
except ImportError:
    import pickle
from collections import OrderedDict
import itertools
import gzip
import numpy
from .utility_functions import timeit


class MLST:
    relative_mlst_path = "strain_typing/resources/pubmlst/mlst_profiles.pkl.gz"

    def __init__(self, plugin_path):
        self.mlst_profile_path = os.path.join(plugin_path, self.relative_mlst_path)
        self.mlst_profiles = pickle.load(gzip.open(self.mlst_profile_path))

    def __repr__(self):
        return str(self.mlst_profiles)

    @timeit
    def find_mlst_profiles(self, strain):
        """


        :param strain:
        :return: found_profiles, status
        """
        results = []
        if self.mlst_profiles is None:
            return None, False

        matching_profiles = OrderedDict()
        matching_counts = OrderedDict()
        matching_sequence = OrderedDict()
        for species, _d in self.mlst_profiles.items():
            matching_profiles.update({species: OrderedDict()})
            matching_counts.update({species: OrderedDict()})
            matching_sequence.update({species: OrderedDict()})
            for i, gene in enumerate(_d["GENE_ORDER"]):
                matching_profiles[species].update({gene: []})
                matching_counts[species].update({gene: {}})
                matching_sequence[species].update({gene: {}})

                for profile_number, profile in _d["GENES"][gene].items():
                    counts = []
                    for kmer in profile:
                        if kmer not in strain.kmers:
                            break
                        if strain.kmers[kmer] < strain.histo_valley_freq:
                            break
                        counts.append(strain.kmers[kmer])
                    else:
                        matching_profiles[species][gene].append(profile_number)
                        matching_counts[species][gene].update({profile_number: counts})
                        matching_sequence[species][gene].update(
                            {profile_number: _d["GENE_SEQUENCE"][gene][profile_number]}
                        )

            st_keys = [":".join(t) for t in list(itertools.product(*matching_profiles[species].values()))]

            for k in st_keys:
                if k in _d["ST"]:
                    if 'Acinetobacter baumannii' in species:
                        st = "{0}_{1}".format(species.split(" ")[-1], _d["ST"][k])
                    else:
                        st = _d["ST"][k]
                else:
                    if 'Acinetobacter baumannii' in species:
                        st = "{0}_{1}".format(species.split(" ")[-1], "NA")
                    else:
                        st = "NA"
                mlst = MLSTProfile(species, st, k.split(":"), _d["GENE_ORDER"])
                mlst.add_counts(matching_counts, matching_sequence)
                results.append(mlst)

        if len(results) == 0:
            return None, True
        return results, True


class MLSTProfile:
    def __init__(self, species, strain_type, profile, gene_order):
        self.species = species
        self.strain_type = strain_type
        self.profile = profile
        self.gene_order = gene_order
        self.allele_counts = []
        self.average_coverage = None
        self.minimum_coverage = None
        self.maximum_coverage = None
        self.sequences = {}

    def print_profile(self):
        s = "species: {0}\n".format(self.species)
        s += "ST: {0}\n".format(self.strain_type)
        s += "profile: {0}\n".format(":".join([str(s) for s in self.profile]))
        s += "genes: {0}\n".format(":".join([str(s) for s in self.gene_order]))
        s += "average_coverage: {0:.2f}\n".format(self.average_coverage) if self.average_coverage is not None else "NA"
        s += "minimum_coverage: {0}\n".format(self.minimum_coverage)
        s += "maximum_coverage: {0}\n".format(self.maximum_coverage)
        for gene, allele_number, kmer_freq in self.allele_counts:
            s += ">{0}-{1}\n{2}\n".format(gene, allele_number, self.sequences[gene])
        for gene, allele_number, kmer_freq in self.allele_counts:
            s += ">{0}-{1}\n{2}\n".format(gene, allele_number, "\t".join(["{0:3}".format(c) for c in kmer_freq]))
        return s

    def add_counts(self, allele_counts, sequences):
        for gene, allele_number in zip(self.gene_order, self.profile):
            self.allele_counts.append((gene, allele_number, allele_counts[self.species][gene][allele_number],))
            self.sequences[gene] = sequences[self.species][gene][allele_number]
        self.average_coverage, self.minimum_coverage, self.maximum_coverage = self.calculate_coverage()

    def calculate_coverage(self):
        all_counts = []
        for gene, allele_number, kmer_freq in self.allele_counts:
            all_counts += kmer_freq
        return numpy.average(all_counts), numpy.min(all_counts), numpy.max(all_counts)

    @staticmethod
    def split_species(species_name):
        try:
            split = species_name.split(" ")
            return "{0}. {1}".format(split[0][0], " ".join(split[1:]))
        except ValueError:
            return "None"

    def __repr__(self):
        if self.average_coverage is None:
            avg_cov = 0
        else:
            avg_cov = int(self.average_coverage)
        return "{0}\tST:{1}\t".format(self.split_species(self.species), self.strain_type,
                                      avg_cov, self.minimum_coverage)


# def main():
#
#     class Strain:
#         def __init__(self, kmers):
#             self.kmers = kmers
#             self.histo_valley_freq = 4
#
#     from os.path import expanduser as eu
#     mlst = MLST("/Users/331-SimmonkLPTP/git_repos/StrainConstructMer")
#     # print(mlst)
#     # load test kmers file 118 (cutoff 45)
#
#     test_file = eu("~/git_repos/StrainConstructMer/strain_typing/resources/test_files/IonCode_0120.kmers.txt")
#     kmer_dict = {}
#     for line in open(test_file):
#         line = line.strip().split("\t")
#         kmer, count = line[0], int(line[1])
#         if count <= 4:
#             continue
#         kmer_dict.update({kmer: count})
#
#     s = Strain(kmer_dict)
#     r = mlst.find_mlst_profiles(s)
#     print(r)
#
#
# if __name__ == "__main__":
#     main()
