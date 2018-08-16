#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os

relation_file = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/get_target.txt"
bed_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/raw_data"
os.chdir(bed_dir)
result_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/02_combined_tf"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


rel_dict = {}
with open(relation_file, 'r') as f:
    for line in f.readlines():
        info = line.strip().split("\t")
        tf = info[0].replace("FLAG-", "").replace("eGFP-", "").replace("-human", "")
        bed_file = info[2].split("/")[-1].split(".gz")[0]
        rel_dict[tf] = rel_dict.get(tf, []) + [bed_file]


def multi_intersect(tissue_name, bed_list):
    in_list = " ".join(bed_list)
    res_file = os.path.join(result_dir, "%s.bed" % tissue_name)
    os.system("cat %s | sort -k1,1 -k2,2n | mergeBed -i - > %s" % (in_list, res_file))
    # os.system("rm %s_*.bed" % tissue_name)


for name, i_list in rel_dict.items():
    if len(i_list) == 1:
        result_file = os.path.join(result_dir, "%s.bed" % name)
        os.system("sort -k1,1 -k2,2n %s | mergeBed -i - > %s" % (i_list[0], result_file))
    else:
        multi_intersect(name, i_list)
