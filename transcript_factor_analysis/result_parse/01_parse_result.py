#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import glob
import collections
from multiprocessing import Pool

result_dir = "/data5/galaxy/project/tf_analysis/interacted_TF/fimo_result"


def get_map_dict():
    total_fimo_files = glob.glob("%s/*/fimo.txt" % result_dir)
    ip_list, input_list, map_dict = [], [], {}
    for fimo_file in total_fimo_files:
        sample_name = fimo_file.split("/")[-2]
        if sample_name.startswith("ip"):
            ip_list.append(fimo_file)
        else:
            input_list.append(fimo_file)
    print(len(ip_list))
    print(len(input_list))
    for ip_file in ip_list:
        tissue_name = ip_file.split("/")[-2].split("_")[1].lower()
        for input_file in input_list:
            i_name = input_file.split("/")[-2].split("_")[0].lower()
            if i_name == tissue_name:
                map_dict[ip_file] = map_dict.get(ip_file, []) + [input_file]
    return map_dict


def combine_all_inputs(input_list):
    all_input_list = []
    for input_file in input_list:
        with open(input_file, 'r') as f:
            f.readline()
            contents = f.readlines()
            contents = [line for line in contents if float(line.split("\t")[-3]) < 0.0001]
            all_input_list += contents
    return all_input_list


def process_single_file(in_list):
    motif_list = [line.split("\t")[0] for line in in_list]
    statistic_result = collections.Counter(motif_list)
    return statistic_result


def calculate_enrich_ratio(ip_file, input_list):
    result_dict = {}
    with open(ip_file, 'r') as f:
        f.readline()
        contents = f.readlines()
        contents = [line for line in contents if float(line.split("\t")[-3]) < 0.0001]
    ip_counter = process_single_file(contents)
    all_input_list = combine_all_inputs(input_list)
    #
    input_counter = process_single_file(all_input_list)
    input_dict = {}
    for i in input_counter:
        input_dict[i] = float(input_counter[i]) / float(len(input_list))
    for i in ip_counter:
        if i in input_dict:
            # correction avoid appear zero
            result_dict[i] = (ip_counter[i] + 0.05) / (input_dict[i] + 0.05)
    #
    tissue_name = ip_file.split("/")[-2].split("_")[1].lower()
    out_file = "%s/%s_fimo.txt" % (result_dir, tissue_name)
    print(out_file)
    with open(out_file, 'w') as fw:
        for i, j in result_dict.items():
            fw.write("%s\t%s\n" % (i, j))


def main_method():
    m_dict = get_map_dict()
    pool = Pool()
    for ip_f, control_list in m_dict.items():
        pool.apply_async(calculate_enrich_ratio, (ip_f, control_list))
    pool.close()
    pool.join()


if __name__ == "__main__":
    main_method()
"""
    for ip_f, inp_list in m_dict.items():
        calculate_enrich_ratio(ip_f, inp_list)
"""