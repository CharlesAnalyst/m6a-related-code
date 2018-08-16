#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import random
import itertools
import pandas as pd
from multiprocessing import Pool


ssap_file = "/data5/galaxy/project/tf_analysis/class_tf/ssap_result/ssap_result.txt"
result_dir = "/data5/galaxy/project/tf_analysis/class_tf/ssap_result"
result_file = os.path.join(result_dir, "parse_result.txt")

# df = pd.read_table(ssap_file, sep="\t", header=None)
# df_sig = df.loc[:, df.loc[:, 4] >= 70]
# print(df_sig.head())


def read_file():
    with open(ssap_file, 'r') as f:
        contents = f.readlines()
    return contents


def get_total_protein(contents):
    pdb_1_list, pdb_2_list = [], []
    for line in contents:
        pdb_1_list.append(line.split()[0]), pdb_2_list.append(line.split()[1])
    total_protein = list(set(pdb_1_list + pdb_2_list))
    return total_protein


def get_relation_dict(contents, pdb_name):
    pdb_dict = {}
    for line in contents:
        if pdb_name in line:
            pdb_1, pdb_2, ssap_score = line.split()[0], line.split()[1], line.split()[2]
            if pdb_1 == pdb_name:
                pdb_dict[pdb_2] = ssap_score
            else:
                pdb_dict[pdb_1] = ssap_score
    return pdb_dict


def get_total_relation_dict(contents, total_protein):
    pdb_dict = {}
    for protein in total_protein:
        i_dict = get_relation_dict(contents, protein)
        pdb_dict[protein] = i_dict
    return pdb_dict


# each pathway represent one cluster;
# one pdb may be clustered into different groups;
def search_pathway(pdb_name, ):
    pass


def main_method():
    contents = read_file()
    total_protein = get_total_protein(contents)
    pdb_dict = get_total_relation_dict(contents, total_protein)
    pass