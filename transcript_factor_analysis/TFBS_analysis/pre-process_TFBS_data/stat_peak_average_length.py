#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd

# 506.254610346bp
####################################################################
m6a_dir = "/data5/galaxy/project/tf_analysis/TFBS_analysis/ip_merge"
m6a_list = glob.glob("%s/*.bed" % m6a_dir)


length_list = []
for m6a in m6a_list:
    print(os.path.basename(m6a))
    df = pd.read_table(m6a, sep="\t", header=None, names=["chr", "start", "end"])
    df["length"] = df["end"] - df["start"]
    length_list.append(df["length"].mean())

print(length_list)
print(sum(length_list) / len(length_list))