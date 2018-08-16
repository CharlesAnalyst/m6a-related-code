#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob


raw_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_hg19_to_GRCh38"
result_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def pick_data():
    index_list = [("chr%d" % i) for i in range(1, 23)]
    index_list.append("chrX")
    index_list.append("chrY")
    for raw_file in glob.glob("%s/*.bed" % raw_dir):
        result_file = os.path.join(result_dir, os.path.basename(raw_file))
        with open(result_file, 'w') as fw:
            with open(raw_file, 'r') as f:
                for line in f.readlines():
                    if line.split("\t")[0] in index_list:
                        fw.write(line)


def sort_data():
    sorted_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean_2"
    for raw_file in glob.glob("%s/*.bed" % result_dir):
        result_file = os.path.join(sorted_dir, os.path.basename(raw_file))
        os.system("sort -k1,1 -k2,2n %s > %s" % (raw_file, result_file))


if __name__ == '__main__':
    pick_data()
    # sort_data()
