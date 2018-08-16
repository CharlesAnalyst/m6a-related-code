#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import time
import pandas as pd


map_file = "/data5/galaxy/project/GWAS_db/00_download_file/gwas_catalog_trait-mappings_r2017-04-24.tsv"
gwas_file = "/data5/galaxy/project/GWAS_db/01_filter_data/filtered_association.txt"
result_dir = "/data5/galaxy/project/GWAS_db/04_class_snp/test"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = "%s/cluster_result.txt" % result_dir


def query_child_parent_relation():
    result_dict, relation_dict, df_gwas = {}, read_map(), read_gwas_file()
    count = 0
    for name, values in df_gwas.iterrows():
        if values["MAPPED_TRAIT"] in relation_dict:
            parent_term = relation_dict[values["MAPPED_TRAIT"]]
            # gwas_bed = "%s\n" % ("\t".join([str(x) for x in values]))
            term = "%s\t%s" % (parent_term, values["MAPPED_TRAIT"])
            # result_dict[parent_term] = result_dict.get(parent_term, []) + [gwas_bed]
            result_dict[term] = result_dict.get(term, 0) + 1
        else:
            print(values["MAPPED_TRAIT"])
    print(count)
    return result_dict


def read_map():
    relation_dict = {}
    df = pd.read_table(map_file, sep="\t")
    trait_list, term_list, parent_term_list = df["Disease trait"], df["EFO term"], df["Parent term"]
    assert len(trait_list) == len(term_list) == len(parent_term_list)
    for i in range(len(trait_list)):
        relation_dict[trait_list[i]] = parent_term_list[i]
        relation_dict[term_list[i]] = parent_term_list[i]
        relation_dict[parent_term_list[i]] = parent_term_list[i]
    return relation_dict


def read_gwas_file():
    df = pd.read_table(gwas_file, sep="\t")
    df["NEW_CHR"], df["START"] = ("chr" + df["CHR_ID"]), (df["CHR_POS"] - 1)
    df_gwas = df[["NEW_CHR", "START", "CHR_POS", "SNPS", "MAPPED_TRAIT"]]
    return df_gwas


def write_to_file(result_dict):
    result_list = []
    for term, number in result_dict.items():
        result_list.append("%s\t%d\n" % (term, number))
    result_list.sort()
    with open(result_file, 'w') as fw:
        fw.writelines(result_list)


if __name__ == "__main__":
    start_time = time.time()
    all_result_dict = query_child_parent_relation()
    write_to_file(all_result_dict)