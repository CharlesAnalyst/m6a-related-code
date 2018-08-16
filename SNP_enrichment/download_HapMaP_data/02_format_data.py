#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import numpy as np
from multiprocessing import Pool
import pandas as pd

#############################################################
vcf_dir = "/data/database/snp"
vcf_list = glob.glob("%s/*.vcf" % vcf_dir)
result_file = "%s/autosomal_GRCh38_snp_unsort.bed" % vcf_dir
final_file = "%s/autosomal_GRCh38_snp_2.bed" % vcf_dir
#############################################################

# memory consume extremely huge!


def process_each_vcf(vcf):
    df = pd.read_table(vcf, sep="\t", comment="#", header=None, names=["CHROM", "POS", "ID", "a", "b", "c", "d", "e"])
    df["chr"], df["start"], df["end"] = ("chr" + df["CHROM"].astype(str)), df["POS"] - 1, df["POS"]
    df_bed = df[["chr", "start", "end", "ID"]].sort_values(by=["chr", "start"]).drop_duplicates()
    print(df_bed.head())
    df_bed.to_csv(result_file, mode="a", sep="\t", header=False, index=False, quoting=False)


def sort_bed():
    os.system("sort -k1,1 -k2,2n %s > %s" % (result_file, final_file))
    # os.remove(result_file)


if __name__ == '__main__':
    """
    pool = Pool(processes=11)
    pool.map_async(process_each_vcf, vcf_list)
    pool.close()
    pool.join()
    #
    print("start sorting")
    sort_bed()
"""
    chr_list = [x.split("chr")[1].split("_")[0] for x in vcf_list]
    index_list = np.argsort(chr_list)
    vcf_list = [vcf_list[i] for i in index_list]
    for i_vcf in vcf_list:
        print(i_vcf)
        process_each_vcf(i_vcf)