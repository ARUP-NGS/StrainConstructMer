from Bio import SeqIO
import os
from .utility_functions import get_kmers
from .utility_functions import timeit
import numpy


class AMRRecord(object):
    def __init__(self, allele, gene_symbol, cds_start, cds_stop, strand, curated_start, gb_nucl_accession,
                 gb_prot_accession, refseq_nucl_accession, proj_accession, description):
        self.allele = allele
        self.gene = gene_symbol
        self.cds_start = cds_start
        self.cds_stop = cds_stop
        self.strand = strand
        self.curated_start = curated_start
        self.gb_nucl_accession = gb_nucl_accession
        self.gb_prot_accession = gb_prot_accession
        self.refseq_nucl_accession = refseq_nucl_accession
        self.proj_accession = proj_accession
        self.description = description.strip()
        self._sequence = None
        self._kmers = None
        self._genus = None
        self._species = None

    def __repr__(self):
        return "{0}\t{1}\t{2}".format(self.allele, self.gene, self.description)

    def hit_formatted(self, kmer_identity, kmer_cutoff, coverage, coverage_change):
        return {
            "allele": self.allele,
            "genus": self.genus,
            "species": self.species,
            "gene": self.gene,
            "description": self.description,
            "accession": self.refseq_nucl_accession,
            "kmer_identity": round(kmer_identity * 100.0, 1),
            "reporting_cutoff": kmer_cutoff,
            "coverage": coverage,
            "coverage_change": coverage_change
        }

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, value):
        self._species = value

    @property
    def genus(self):
        return self._genus

    @genus.setter
    def genus(self, value):
        self._genus = value

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        self._sequence = value.upper()

    @property
    def kmers(self):
        return self._kmers

    @kmers.setter
    def kmers(self, value):
        self._kmers = set(value)


class AMRDatabase:
    amr_meta_path = 'strain_typing/resources/AMR/Bacterial_Antimicrobial_Resistance_Reference_Gene_Database_020718a.txt'
    arm_sequence_path = 'strain_typing/resources/AMR/amr.fasta'

    def __init__(self, plugin_dir):
        self.amr_records = self.__parse_amr_descriptions(plugin_dir)
        self.add_sequence_to_record(plugin_dir)

    @timeit
    def add_sequence_to_record(self, plugin_dir):
        for s in SeqIO.parse(os.path.join(plugin_dir, self.arm_sequence_path), 'fasta'):
            accession = s.name
            genus = s.description.split(" ")[1]
            species = " ".join(s.description.split(" ")[1:3])
            if accession not in self.amr_records:
                print(accession)
            self.amr_records[accession].sequence = str(s.seq)
            self.amr_records[accession].kmers = get_kmers(str(s.seq))
            self.amr_records[accession].genus = genus
            self.amr_records[accession].species = species

    @timeit
    def __parse_amr_descriptions(self, plugin_dir):
        records = {}
        for line in open(os.path.join(plugin_dir, self.amr_meta_path)):
            if line[0] == "#":
                continue

            line = line.split("\t")
            allele, gene_symbol, cds_start, cds_stop, strand, curated_start, gb_nucl_accession, gb_prot_accession, \
            refseq_nucl_accession, proj_accession, description = line

            if allele == "":
                allele = gene_symbol

            if gene_symbol == "":
                gs = description.split("\t")[-1]
                allele = gs
                gene_symbol = gs

            arm_record = AMRRecord(allele, gene_symbol, cds_start, cds_stop, strand, curated_start, gb_nucl_accession,
                                   gb_prot_accession, refseq_nucl_accession, proj_accession, description)
            records[refseq_nucl_accession] = arm_record
        return records

    @staticmethod
    def __get_kmer_identity(gene_kmers, strain_kmers):
        matched = 0
        total = 0
        coverage = []
        for kmer in gene_kmers:
            total += 1
            if kmer in strain_kmers:
                matched += 1
                coverage.append(strain_kmers[kmer])
            else:
                coverage.append(0)
        return float("{0:.3f}".format(float(matched)/total)), coverage

    @timeit
    def find_antibiotic_genes(self, strain, coverage_cutoff=.80):
        amr_hits = []
        for amr_record in self.amr_records.values():
            kid, coverage = self.__get_kmer_identity(amr_record.kmers, strain.kmers)
            # FAILED
            if kid < coverage_cutoff:
                continue
            # PASSES
            diff_cov = self.calculate_coverage_difference(strain.coverage, coverage)
            formatted_record = amr_record.hit_formatted(kid, coverage_cutoff, coverage, diff_cov)
            amr_hits.append(formatted_record)

        amr_hits = sorted(amr_hits, key=lambda k: k["kmer_identity"], reverse=True)

        return amr_hits  # this contains duplicates genes which are reduced in the strain.json output

    @staticmethod
    def calculate_coverage_difference(strain_coverage, gene_coverage):
        # return float("{0:.1f}".format(float(numpy.mean(numpy.nonzero(gene_coverage))) / strain_coverage))
        return float("{0:.1f}".format(float(numpy.median([i for i in gene_coverage if i > 0])) / strain_coverage))
