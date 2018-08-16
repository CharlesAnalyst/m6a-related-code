#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import pandas as pd

# select transcript factor must be in enrich tf(m6a)
###################################################################################################
tf_file = "/data5/galaxy/project/tf_pathway/mettl3_14_watp/TF_symbol.txt"
interaction_file_1 = "/data4/database/gene_interaction/export_Tue_Feb_06_07_42_10_UTC_2018.csv"
interaction_file_2 = "/data5/galaxy/project/tf_pathway/mettl3_14_watp/methyl_interacted_protein.txt"
result_dir = "/data5/galaxy/project/tf_pathway/mettl3_14_watp"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
############################################


df = pd.read_table(tf_file, sep="\t", header=None, names=["transcript factor"])
tf_list = df["transcript factor"].tolist()

result_file, result_list = "%s/all_interaction.txt" % result_dir, []
with open(interaction_file_1, 'r') as f:
    f.readline()
    contents = f.readlines()
    for line in contents:
        info = line.strip().split(",")
        if info[0] in tf_list and info[2] in tf_list:
            result_line = "%s\tinteraction\t%s" % (info[0], info[2])
            print(result_line)
            result_list.append(result_line)
with open(interaction_file_2, 'r') as f:
    f.readline()
    for line in contents:
        info = line.strip().split("\t")
        if info[0] in tf_list and info[2] in tf_list:
            result_line = line.strip()
            result_list.append(result_line)

uniq_list = list(set(result_list))
with open(result_file, 'w') as fw:
    fw.write("source\tact\ttarget\n")
    for line in uniq_list:
        fw.write(line + "\n")
