#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import collections
import pandas as pd
from multiprocessing import Pool

############################################################################
result_dir = "/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result"
tmp_dir = "%s/combined_input" % result_dir
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)
input_cycle_num = 100.0
############################################################################


def get_ip_and_input_files():
    total_fimo_files = glob.glob("%s/*/fimo.txt" % result_dir)
    ip_dict, input_dict = {}, {}
    for fimo_file in total_fimo_files:
        sample_name = fimo_file.split("/")[-2]
        if sample_name.startswith("ip"):
            tissue_name = sample_name.split("_")[1].lower()
            ip_dict[tissue_name] = fimo_file
        else:
            tissue_name = sample_name.split("_")[0].lower()
            input_dict[tissue_name] = input_dict.get(tissue_name, []) + [fimo_file]
    return ip_dict, input_dict


def get_combined_input_dict(input_dict):
    result_dict = {}
    for tissue, input_list in input_dict.items():
        print(tissue, input_list[0])
        combined_file = os.path.join(tmp_dir, "%s_fimo.txt" % input_list[0].split("/")[-2].split("_")[0].lower())
        # os.system("awk 'FNR>1' %s | awk '{print $1}' > %s" % (" ".join(input_list), combined_file))
        result_dict[tissue] = combined_file
    return result_dict


def process_single_file(in_file):
    df = pd.read_table(in_file, sep="\t", comment="#", header=None)
    motif_list = df.loc[:, 0].tolist()
    statistic_result = collections.Counter(motif_list)
    return statistic_result


def calculate_enrich_ratio(ip_file, input_file):
    result_dict, input_dict = {}, {}
    ip_counter = process_single_file(ip_file)
    input_counter = process_single_file(input_file)
    for i in input_counter:
        input_dict[i] = float(input_counter[i]) / input_cycle_num
    for i in ip_counter:
        if i in input_dict:
            result_dict[i] = (ip_counter[i] + 0.05) / (input_dict[i] + 0.05)
    tissue_name = ip_file.split("/")[-2].split("_")[1].lower()
    out_file = "%s/%s_fimo.txt" % (result_dir, tissue_name)
    print(out_file)
    with open(out_file, 'w') as fw:
        for i, j in result_dict.items():
            fw.write("%s\t%s\n" % (i, j))


def main_method():
    ip_dict, input_dict = get_ip_and_input_files()
    input_dict = get_combined_input_dict(input_dict)
    pool = Pool()
    for tissue, ip_file in ip_dict.items():
        input_file = input_dict[tissue]
        print(ip_file, input_file)
        pool.apply_async(calculate_enrich_ratio, (ip_file, input_file))
    pool.close()
    pool.join()
    # os.system("rm -rf %s" % tmp_dir)


if __name__ == "__main__":
    main_method()
