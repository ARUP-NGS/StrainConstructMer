import sys
import gzip
try:
    import cPickle as pickle
except ImportError:
    import pickle

from Bio import SeqIO
import os

from collections import OrderedDict
import requests


mlst_urls = {
    'Acinetobacter baumannii Oxf': [
        'http://pubmlst.org/data/profiles/abaumannii.txt',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_gltA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_gyrB.tfa',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_gdhB.tfa',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_recA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_cpn60.tfa',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_gpi.tfa',
        'http://pubmlst.org/data/alleles/abaumannii/Oxf_rpoD.tfa',
    ],
    'Acinetobacter baumannii Pas': [
        'http://pubmlst.org/data/profiles/abaumannii_2.txt',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_cpn60.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_fusA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_gltA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_pyrG.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_recA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_rplB.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_rpoB.tfa',
    ],
    'Enterococcus faecalis': [
        'http://pubmlst.org/data/profiles/efaecalis.txt',
        'http://pubmlst.org/data/alleles/efaecalis/gdh.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/gyd.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/pstS.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/gki.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/aroE.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/xpt.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/yqiL.tfa',
    ],
    'Enterococcus faecium': [
        'http://pubmlst.org/data/profiles/efaecium.txt',
        'http://pubmlst.org/data/alleles/efaecium/atpA.tfa',
        'http://pubmlst.org/data/alleles/efaecium/ddl.tfa',
        'http://pubmlst.org/data/alleles/efaecium/gdh.tfa',
        'http://pubmlst.org/data/alleles/efaecium/purK.tfa',
        'http://pubmlst.org/data/alleles/efaecium/gyd.tfa',
        'http://pubmlst.org/data/alleles/efaecium/pstS.tfa',
        'http://pubmlst.org/data/alleles/efaecium/adk.tfa',
    ],
    'Staphylococcus aureus': [
        'http://pubmlst.org/data/profiles/saureus.txt',
        'http://pubmlst.org/data/alleles/saureus/arcC.tfa',
        'http://pubmlst.org/data/alleles/saureus/aroE.tfa',
        'http://pubmlst.org/data/alleles/saureus/glpF.tfa',
        'http://pubmlst.org/data/alleles/saureus/gmk.tfa',
        'http://pubmlst.org/data/alleles/saureus/pta_.tfa',
        'http://pubmlst.org/data/alleles/saureus/tpi.tfa',
        'http://pubmlst.org/data/alleles/saureus/yqiL.tfa',
    ],
    'Salmonella enterica': [
        'http://pubmlst.org/data/profiles/senterica.txt',
        'http://pubmlst.org/data/alleles/senterica/aroC.tfa',
        'http://pubmlst.org/data/alleles/senterica/dnaN.tfa',
        'http://pubmlst.org/data/alleles/senterica/hemD.tfa',
        'http://pubmlst.org/data/alleles/senterica/hisD.tfa',
        'http://pubmlst.org/data/alleles/senterica/purE.tfa',
        'http://pubmlst.org/data/alleles/senterica/sucA.tfa',
        'http://pubmlst.org/data/alleles/senterica/thrA.tfa',
    ],
    'Escherichia coli': [
        'http://pubmlst.org/data/profiles/ecoli.txt',
        'http://pubmlst.org/data/alleles/ecoli/adk.tfa',
        'http://pubmlst.org/data/alleles/ecoli/fumC.tfa',
        'http://pubmlst.org/data/alleles/ecoli/gyrB.tfa',
        'http://pubmlst.org/data/alleles/ecoli/icd.tfa',
        'http://pubmlst.org/data/alleles/ecoli/mdh.tfa',
        'http://pubmlst.org/data/alleles/ecoli/purA.tfa',
        'http://pubmlst.org/data/alleles/ecoli/recA.tfa'
    ]
}


def reverse_complement(sequence):
    """
    return the reverse complement of the sequence string
    :param sequence:
    :return:
    """
    base_map = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "Y": "R", "R": "Y", "M": "K", "K": "M",
                "S": "W", "W": "S", "V": "D", "D": "V", "H": "B", "B": "H"}
    try:
        new_sequence = "".join([base_map[char] for char in reversed(sequence)])
    except KeyError as e:
        raise KeyError("{0} char not present in sequence".format(e))
    return new_sequence


def update_mlst_resources():
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    # this is the path to save resources
    resource_path = os.path.join(_ROOT, '..', 'strain_typing', 'resources', 'pubmlst')
    if os.path.exists(resource_path) is False:
        os.mkdir(resource_path)
    file_lists = {}
    for k in mlst_urls.keys():
        strain_path = os.path.join(resource_path, k)
        if os.path.exists(strain_path) is False:
            os.mkdir(strain_path)

        sys.stderr.write("Updating '{0}' resources\n".format(k))
        sys.stderr.write("...saving_to '{0}' resources\n".format(resource_path))
        file_lists.update({k: []})
        for url in mlst_urls[k]:
            file_lists[k].append(os.path.join(strain_path, url.split("/")[-1].lower().replace("_.", ".")))
            sys.stderr.write('\tretrieving: {0}\n'.format(url))
            sys.stderr.write('\tsaving: {0}\n'.format(url.split("/")[-1].lower().replace("_.", ".")))
            with open(os.path.join(strain_path, url.split("/")[-1].lower().replace("_.", ".")), 'wb') as wf:
                content = requests.get(url, verify=False)
                wf.write(content.content)


    pickle_profiles(file_lists, resource_path)
    return


def pickle_profiles(file_lists, resource_path, kmer_size=31):
    # jellyfish.MerDNA_k(KMER_SIZE)
    # instantiate the pickle obj
    mlst_profiles_dict = OrderedDict()
    for species, file_list in file_lists.items():
        if species not in mlst_profiles_dict:
            mlst_profiles_dict.update({species: {"ST": OrderedDict(), "GENES": OrderedDict(),
                                                 "GENE_SEQUENCE": OrderedDict(), "GENE_ORDER": None}})

        number_of_genes = len(file_list) - 1
        for i, l in enumerate(open([os.path.join(species, f)
                                    for f in file_list if '.txt' in f][0])):
            line = l.strip().split("\t")
            if i == 0:
                gene_list = line[1:number_of_genes + 1]
            else:
                profile = ":".join(line[1:number_of_genes + 1])
                st = line[0]
                mlst_profiles_dict[species]["ST"].update({profile: st})
                mlst_profiles_dict[species]["GENE_ORDER"] = gene_list

        for _file in [f for f in file_list if f[-4:] == '.tfa']:
            for seq_record in SeqIO.parse(os.path.join(resource_path, species, _file), 'fasta'):
                seq_num = seq_record.name.replace("-", "_").split("_")[-1]
                gene_name = "_".join(seq_record.name.replace("__", "_").replace("-", "_").split("_")[:-1])

                if gene_name not in mlst_profiles_dict[species]["GENES"]:
                    mlst_profiles_dict[species]["GENES"].update({gene_name: {seq_num: set([])}})
                else:
                    mlst_profiles_dict[species]["GENES"][gene_name].update({seq_num: set([])})

                if gene_name not in mlst_profiles_dict[species]["GENE_SEQUENCE"]:
                    mlst_profiles_dict[species]["GENE_SEQUENCE"].update({gene_name: {seq_num: str(seq_record.seq)}})

                else:
                    mlst_profiles_dict[species]["GENE_SEQUENCE"][gene_name].update({seq_num: str(seq_record.seq)})

                for j in range(0, len(seq_record.seq) - kmer_size + 1):
                    kmer = str(seq_record.seq[j: j + kmer_size])
                    mlst_profiles_dict[species]["GENES"][gene_name][seq_num].add(min(kmer, reverse_complement(kmer)))
            sys.stderr.write("\tparsing: {0} : {1}\n".format(species, gene_name))

    with gzip.GzipFile(os.path.join(resource_path, 'mlst_profiles.pkl.gz'), "w") as wf:
        pickle.dump(mlst_profiles_dict, wf)
    return


def main():
    update_mlst_resources()


if __name__ == '__main__':
    main()
