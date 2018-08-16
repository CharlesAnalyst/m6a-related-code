#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob


data_dir = ""


os.chdir("/data5/galaxy/project/trend_test")


def add_group_info(in_file):
    prefix = in_file.split("_")[0]
    rela_dict = {"%s-1" % prefix: "low", "%s-2" % prefix: "low", "%s-3" % prefix: "medium", "%s-4" % prefix: "medium", "%s-5" % prefix: "medium", "%s-6" % prefix: "high", "%s-7" % prefix: "high", "%s-8" % prefix: "high"}
    result_file = "group_%s" % in_file
    with open(result_file, 'w') as fw:
        with open(in_file, 'r') as f:
            for line in f.readlines():
                a, b = line.strip().split("\t")
                new_line = "%s\t%s\t%s\t%s\n" % (a, b, rela_dict[b], prefix)
                fw.writelines(new_line)


# for file in ["total_data.txt", "mRNA_data.txt", "linc_data.txt"]:
for file in ["5utr_data.txt", "cds_data.txt", "stop_data.txt", "3utr_data.txt"]:
    add_group_info(file)