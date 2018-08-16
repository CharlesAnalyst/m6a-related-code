#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd

##############################################################################################
tf_dir = "/data5/galaxy/project/tf_pathway/m6aSNP_TF_intersect/TF"
# tf_file_list = glob.glob("%s/*.txt" % tf_dir)
# export csv from DisGeNet
interaction_file = "/data4/database/gene_interaction/export_Tue_Feb_06_07_42_10_UTC_2018.csv"
result_dir = "/data5/galaxy/project/tf_pathway/gene_interaction"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
############################################


os.chdir(tf_dir)
df = pd.read_table("all_tissues.txt", sep="\t")
tf_list = df["transcript factor"].tolist()
# tf_list = [tf_list[i] for i in range(10)]

result_file = "%s/all_tissues.txt" % result_dir
with open(result_file, 'w') as fw:
    fw.write("source\tact\ttarget\n")
    with open(interaction_file, 'r') as f:
        f.readline()
        contents = f.readlines()
        for line in contents:
            info = line.strip().split(",")
            if info[0] in tf_list and info[2] in tf_list:       # select!
                fw.write("%s\tinteraction\t%s\n" % (info[0], info[2]))
