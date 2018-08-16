#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import sys
import glob
import pandas as pd

data_file = "/data5/galaxy/project/tf_analysis/class_tf/CATH_db/human_annotations.tsv"
cath_file = "/data5/galaxy/project/tf_analysis/class_tf/CATH_db/cath-names.txt"
tfbs_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/02_combined_tf"
result_file = "/data5/galaxy/project/tf_analysis/class_tf/CATH_db/cluster_result_function.txt"


def class_cf():
    name_dict, pdb_dict, gene_list = get_cath_names_dict(), get_super_family(), get_TFBS_gene_names()
    fw = open(result_file, 'w')
    fw.write("PDB Code\tDescription\tGene number\tTFBS number\tTFBS name\n")
    for pdb_code, genes in pdb_dict.items():
        genes = list(set(genes))
        try:
            description = name_dict[pdb_code]
            tfbs_gene = list(set([x for x in genes if x in gene_list]))
            fw.write("%s\t%s\t%d\t%d\t%s\n" % (pdb_code, description, len(genes), len(tfbs_gene), ";".join(tfbs_gene)))
        except KeyError:
            print("%s not in name_dict!" % pdb_code)
    fw.close()


def get_super_family():
    pdb_dict = {}
    df = pd.read_table(data_file, sep="\t")
    df_human = df[df["TAXON_ID"] == 9606]
    for name, values in df_human.iterrows():
        i_key = values["SUPERFAMILY"]
        pdb_dict[i_key] = pdb_dict.get(i_key, []) + [values["GENE_ID"].split("_HUMAN")[0]]
    return pdb_dict


def get_cath_names_dict():
    name_dict = {}
    df = pd.read_table(cath_file, sep="\s\s\s\s", comment="#", header=None, names=["SUPERFAMILY", "pdb", "description"])
    print(df.head())
    for name, values in df.iterrows():
        name_dict[values["SUPERFAMILY"]] = values["description"].replace(":", "")
    return name_dict


def get_TFBS_gene_names():
    gene_list = [os.path.basename(x).split(".bed")[0].upper() for x in glob.glob("%s/*.bed" % tfbs_dir)]
    return gene_list


if __name__ == '__main__':
    class_cf()
