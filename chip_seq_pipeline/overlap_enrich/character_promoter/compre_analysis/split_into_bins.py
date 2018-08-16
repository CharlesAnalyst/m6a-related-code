#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import pandas as pd


############################
UPSTREAM = 3000
DOWNSTREAM = 3000
PROMOTER_UPSTREAM = 2000
BIN_SIZE = 200
#############
promoter_bed = "/data5/galaxy/project/data/promoter/human/2k_100/genes_promoters.bed"
tss_bed = "/data5/galaxy/project/mettl3_enrich/mettl3_enrich_high_CpG/tss_high.bed"
#############################


def main():
    # df_tss_region = get_tss_region()
    df_tss_region = pd.read_table(tss_bed, sep="\t", header=None, names=["chr", "tss_s", "tss_e", "d", "e", "strand"])
    bin_list = split_tss_region_into_bins(df_tss_region)
    return bin_list


def get_tss_region():
    df = pd.read_table(promoter_bed, sep="\t", header=None, names=["chr", "start", "end", "d", "e", "strand"])
    df_pos, df_neg = df[df.loc[:, "strand"] == "+"], df[df.loc[:, "strand"] == "-"]
    #
    df_neg["tss_s"], df_neg["tss_e"] = ((df_neg["end"] - PROMOTER_UPSTREAM) - DOWNSTREAM), (
            (df_neg["end"] - PROMOTER_UPSTREAM) + UPSTREAM)
    df_neg_bed = df_neg[["chr", "tss_s", "tss_e", "d", "e", "strand"]]
    #
    df_pos["tss_s"], df_pos["tss_e"] = ((df_pos["start"] + PROMOTER_UPSTREAM) - UPSTREAM), (
            (df_pos["start"] + PROMOTER_UPSTREAM) + DOWNSTREAM)
    df_pos_bed = df_pos[["chr", "tss_s", "tss_e", "d", "e", "strand"]]
    #
    df_tss_region = pd.concat([df_pos_bed, df_neg_bed]).sort_values(["chr", "tss_s"])
    df_tss_region = df_tss_region[df_tss_region["tss_s"] > 0]
    return df_tss_region


def split_tss_region_into_bins(df_tss_region):
    bin_list = []   # Did not keep strand information!
    for name, values in df_tss_region.iterrows():
        chro, start, end, each_bin_list = values["chr"], int(values["tss_s"]), int(values["tss_e"]), []
        for i in range(start, end, BIN_SIZE):
            each_bin_list.append("%s\t%d\t%d" % (chro, i, i+BIN_SIZE))
        bin_list.append(each_bin_list)
    return bin_list

