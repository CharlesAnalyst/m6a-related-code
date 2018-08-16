#!/usr/bin/python
# -*- coding: utf-8 -*-
# @author galaxy

import os
import sys
import glob
import subprocess
import pandas as pd

#####################################################################################################
snp_bed = "/data5/galaxy/project/snp_analysis/GTEx_analysis/gain_or_loss/02_parse_result/GRCh38.bed"
m6a_dir = "/data5/galaxy/project/data/total_m6a_peak"
m6a_list = glob.glob("%s/*.bed" % m6a_dir)
tissue_list = ["BRAIN", "HEART", "LIVER", "LUNG", "MUSCLE", "STOMACH"]
######
snp_db = "/data/database/snp/snp150Common_ucsc_hg38/snp.bed"
df_snp = pd.read_table(snp_db, sep="\t", header=None, names=["chr", "start", "end", "id", "s", "strand"])
snp_dict = dict(zip(df_snp["end"], df_snp["id"]))
print(snp_dict[225574923])
######
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/create_big-table"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
######################################################################################################


def get_snp_content():
    df = pd.read_table(snp_bed, sep="\t", header=None, names=["chr", "start", "end", "a", "b"])
    df["tissue"] = df["a"].str.split(";").str[0]
    df["function"] = df["b"].str.split(";").str[0]
    df["ref"] = df["b"].str.split("|").str[-3]
    df["alter"] = df["b"].str.split("|").str[-2]
    df["other"] = df["tissue"].astype(str) + ";" + df["function"].astype(str) + ";" + df["ref"].astype(str) + ";" + df["alter"].astype(str)
    df_group = df.groupby(["tissue"]).apply(match_function)


def match_function(df):
    pro_tissue = df.iloc[0, -5]
    tissue = df.iloc[0, -5].split("_")[0].upper()
    print(tissue)
    m6a_dict = get_each_tissue_files(m6a_list)
    m6a_bed = m6a_dict[tissue][0]
    print(m6a_bed)
    result_file = os.path.join(result_dir, "%s.txt" % pro_tissue)
    df_result = intersect_with_m6a(df, m6a_bed)
    # print(df_result.head())
    df = get_rs_id(df_result)
    print(df.head())
    df.to_csv(result_file, sep="\t", header=None, index=False)


def get_each_tissue_files(file_list):
    tissue_file_dict = {}
    for i_file in file_list:
        tissue = os.path.basename(i_file).split(".bed")[0].split("_")[0].upper()
        tissue_file_dict[tissue] = tissue_file_dict.get(tissue, []) + [i_file]
    return tissue_file_dict


def intersect_with_m6a(df, m6a_bed):
    # print(df.head())
    df = df.sort_values(["chr", "start", "end"])
    df_str = df[["chr", "start", "end", "other"]].to_string(header=False, index=False, col_space=4)
    snp_str = "\n".join(["\t".join(x.split()) for x in df_str.split("\n")])
    snp_str = snp_str.encode("utf-8")
    # print(snp_str)
    command = "bedtools intersect -a stdin -b %s -wa -wb" % m6a_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    results = sub_p.communicate(snp_str)[0].decode("utf-8").strip().split("\n")
    # print(results)
    n_arrays = [x.split("\t") for x in results]
    df_result = pd.DataFrame(n_arrays)
    return df_result


def get_rs_id(df):
    df["tissue"] = df.iloc[:, 3].str.split(";").str[0]
    df["function"] = df.iloc[:, 3].str.split(";").str[1]
    df["ref"] = df.iloc[:, 3].str.split(";").str[2]
    df["alter"] = df.iloc[:, 3].str.split(";").str[3]
    df["rs id"] = [snp_dict[int(x)] if (int(x) in snp_dict) else x for x in df.iloc[:, 2]]
    df_final = df[["rs id", "ref", "alter", "function", "tissue"]]
    df_final = df_final.drop_duplicates()
    return df_final


if __name__ == '__main__':
    get_snp_content()