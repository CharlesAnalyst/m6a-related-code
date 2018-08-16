#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

# size match: total length of chromosome
import os
import datetime
import pandas as pd
from multiprocessing import Pool
from sklearn.utils import shuffle


###########################################################################################################
cycle_number = 100
m6a_bed = "/data5/galaxy/project/data/total_m6a_peak/universal/universal.bed"
control_bed = "/data5/galaxy/project/data/input_data/autosomal/universal/input_bed/universal_discrete.bed"
result_dir = "/data5/galaxy/project/GWAS_analysis/binomial_test/gwas_vs_background/Fisher/100cycle_input"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
############################################################################################################


def pick_control_peak():
    chromosome_len_dict = count_ip_length_byChrom()
    df = pd.read_table(control_bed, sep="\t", header=None, comment="#", names=["chr", "start", "end"])
    df.loc[:, "length"] = df["end"] - df["start"]
    pool = Pool()
    for i in range(cycle_number):
        df = shuffle(df)
        pool.apply_async(pick_control_length_byChrom, (i+1, df, chromosome_len_dict))
        # pick_control_length_byChrom(i + 1, df, chromosome_len_dict)
    pool.close()
    pool.join()


def count_ip_length_byChrom():
    chromosome_len_dict, df_list = {}, []
    df = pd.read_table(m6a_bed, sep="\t", header=None, comment="#", names=["chr", "start", "end"])
    df.loc[:, "length"] = df["end"] - df["start"]
    chr_list = list(set(df.loc[:, "chr"].tolist()))
    for chr_name in chr_list:
        df_i = df[df["chr"] == chr_name]
        chr_sum_len = df_i["length"].sum()
        chromosome_len_dict[chr_name] = chr_sum_len
    return chromosome_len_dict


def pick_control_length_byChrom(index, df, chr_len_dict):
    total_list = []
    chr_list = list(chr_len_dict.keys())
    for chr_name in chr_list:
        df_chr = df[df["chr"] == chr_name]
        chr_length = chr_len_dict[chr_name]
        peak_list = pick_equal_length_byChrom(chr_length, df_chr)
        total_list += peak_list
    df_result = pd.DataFrame(total_list, columns=["chr", "start", "end"])
    result_file = os.path.join(result_dir, "control-%d.bed" % index)
    df_sorted = df_result.sort_values(["chr", "start"])
    df_sorted.to_csv(result_file, sep="\t", header=None, index=False)
    

def pick_equal_length_byChrom(chr_length, df_chr):
    df_chr.columns = ["chr", "start", "end", "length"]
    count, peak_list = 0, []
    if df_chr["length"].sum() < chr_length:
        print("Input chr length less than IP!")
    for name, values in df_chr.iterrows():
        peak_list.append([values["chr"], values["start"], values["end"]])
        count += values["length"]
        if count >= chr_length:
            break
    chr_name, start, end = peak_list[-1]
    chr_name, start, new_end = chr_name, int(start), (int(end) - (count - chr_length))
    peak_list[-1] = [chr_name, start, new_end]
    return peak_list


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    pick_control_peak()
    print(start_time)
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

