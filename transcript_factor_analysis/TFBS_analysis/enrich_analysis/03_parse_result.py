#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
import numpy as np
from collections import Counter
from multiprocessing import Pool


####################################################################################
in_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/04_intersect"
tissue_list = ["brain", "heart", "kidney", "liver", "lung", "muscle", "placenta", "stomach"]
# in_file_list = glob.glob("%s/*/fimo.txt" % in_dir)
##########################################################


def statistic_number_each_tissue(tissue_file):
    result_file = tissue_file.replace("fimo.txt", "statistic_result.txt")
    df = pd.read_table(tissue_file, sep="\t")
    # print(df.head())
    peak_list = list(set(df["sequence_name"]))
    total_list = []
    for peak in peak_list:
        df_sub = df[df["sequence_name"] == peak]
        motif_list = df_sub["# motif_id"]
        counter = Counter(motif_list)
        result_list = [[peak, i_motif, counter[i_motif]] for i_motif in counter if counter[i_motif] >= 3]
        total_list += result_list
    print(total_list[0])
    df_result = pd.DataFrame(total_list, columns=["peak", "motif", "counter"])
    print(df_result.head())
    df_sort = df_result.sort_values(["peak", "motif", "counter"])
    df_sort.to_csv(result_file, sep="\t", index=False, quoting=False)


def statistic_file_line(in_file):
    tmp = os.path.basename(in_file).split("-")
    if len(tmp) == 2:
        prefix = tmp[1].split(".bed")[0]
    else:
        prefix = tmp[0].split(".bed")[0]
    return_s = os.popen("wc -l %s" % in_file)
    overlap_num = return_s.read().split()[0]
    return "%s\t%s" % (prefix, str(overlap_num))


def statistic_total_tissue():
    for tissue in tissue_list:
        print(tissue)
        ip_dir, input_dir = os.path.join(in_dir, tissue, "ip"), os.path.join(in_dir, tissue, "control")
        ip_count_dict, input_count_dict = {}, {}
        ip_list, input_list = glob.glob("%s/*" % ip_dir), glob.glob("%s/*" % input_dir)
        for ip in ip_list:
            result = statistic_file_line(ip)
            ip_count_dict[result.split("\t")[0]] = int(result.split("\t")[1])
        pool, result_list = Pool(), []
        for in_file in input_list:
            result = pool.apply_async(statistic_file_line, (in_file, ))
            result_list.append(result)
        pool.close()
        pool.join()
        for result in result_list:
            prefix, num = result.get().split("\t")
            input_count_dict[prefix] = input_count_dict.get(prefix, []) + [int(num)]
        result_file = "%s/%s.txt" % (os.path.join(in_dir, tissue), tissue)
        with open(result_file, 'w') as fw:
            for prefix, ip_num in ip_count_dict.items():
                input_num_list = input_count_dict[prefix]
                input_num_mean = np.mean(input_num_list)
                enrich_ratio = (ip_num + 0.05) / (input_num_mean + 0.05)
                fw.write("%s\t%f\n" % (prefix, enrich_ratio))


if __name__ == "__main__":
    statistic_total_tissue()
