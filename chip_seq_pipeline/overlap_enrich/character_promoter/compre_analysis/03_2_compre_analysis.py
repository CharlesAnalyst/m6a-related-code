#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import subprocess
import pandas as pd
import split_into_bins


#############################################################################################
raw_mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3_peaks.narrowPeak"
mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
result_dir = "/data5/galaxy/project/mettl3_enrich/mettl3_enrich_high_CpG/peak_version"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
# result_file = os.path.join(result_dir, "mettl3_vs_free.txt")
#######################################################################


def main():
    bin_list = split_into_bins.main()
    peak_str_list = get_occupied_promoter(bin_list)
    for peak_str in peak_str_list:
        df = pd.DataFrame([peak.split("\t") for peak in peak_str.strip().split("\n")])
        df.columns = ["chr", "start", "end"]
        df["start"], df["end"] = df["start"].astype(int), df["end"].astype(int)
        average_pvalue = query_peak_pvalue(df)
        print(average_pvalue)


def get_occupied_promoter(bin_list):
    peak_str_list = []
    for i in range(len(bin_list[0])):
        each_pos_bin_list = []
        for j in range(len(bin_list)):
            each_pos_bin_list.append(bin_list[j][i])
        each_pos_bin_str = "\n".join(each_pos_bin_list)
        str_peak = overlap_with_peak_each_bin(each_pos_bin_str)
        peak_str_list.append(str_peak)
    return peak_str_list


def overlap_with_peak_each_bin(each_pos_bin_str):
    command = "bedtools intersect -a %s -b stdin -wa | sort -k1,1 -k2,2n | uniq " % mettl3_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    str_peak = sub_p.communicate(each_pos_bin_str)[0]
    return str_peak


def query_peak_pvalue(df_peak):
    df = pd.read_table(raw_mettl3_bed, sep="\t", header=None, names=["chr", "start", "end", "d", "e", "f", "g", "pvalue", "i", "j"])
    df_com = pd.merge(df_peak, df, on=["chr", "start", "end"], how="left").dropna(how="any")
    df_com = df_com[["chr", "start", "end", "d", "pvalue", "f"]]
    # print(df_com)
    average_pvalue = df_com["pvalue"].sum()
    return average_pvalue


if __name__ == '__main__':
    main()



































