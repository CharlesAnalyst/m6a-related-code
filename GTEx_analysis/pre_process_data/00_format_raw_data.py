#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd

#########################################
file_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/raw_data"
file_list = [file for file in glob.glob("%s/*.txt" % file_dir) if "signif_variant" in file]
#######
vcf_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/01_1_format_vcf_hg19"
bed_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/01_2_format_bed_GRCh37"
for i_dir in [vcf_dir, bed_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)
############################################


def generate_vcf():
    for i_file in file_list:
        prefix = os.path.basename(i_file).split(".")[0]
        result_file = os.path.join(vcf_dir, "%s.txt" % prefix)
        df = pd.read_table(i_file, sep="\t")
        df["sample"] = prefix
        df["chr"] = df["variant_id"].str.split("_").str[0]
        df["pos"] = df["variant_id"].str.split("_").str[1].astype(int)
        df["new_chr"] = "chr" + df["chr"].astype(str)
        df["new_pos"] = df["pos"] - 1
        df["ref"] = df["variant_id"].str.split("_").str[2]
        df["alt"] = df["variant_id"].str.split("_").str[3]
        df_vcf = df[["new_chr", "new_pos", "variant_id", "ref", "alt", "sample"]]
        df_sort = df_vcf.sort_values(["new_chr", "new_pos", "variant_id"])
        df_sort.to_csv(result_file, sep="\t", header=None, index=False)


def generate_bed():
    for i_file in file_list:
        prefix = os.path.basename(i_file).split(".")[0]
        result_file = os.path.join(bed_dir, "%s.bed" % prefix)
        df = pd.read_table(i_file, sep="\t")
        # df["sample"] = prefix
        df["chr"] = df["variant_id"].str.split("_").str[0]
        df["pos"] = df["variant_id"].str.split("_").str[1].astype(int)
        df["new_chr"] = "chr" + df["chr"].astype(str)
        df["new_pos"] = df["pos"] - 1
        df_vcf = df[["new_chr", "new_pos", "pos", "variant_id"]]
        df_sort = df_vcf.sort_values(["new_chr", "new_pos", "variant_id"])
        df_sort.to_csv(result_file, sep="\t", header=None, index=False)


if __name__ == '__main__':
    generate_vcf()
    # generate_bed()
