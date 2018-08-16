#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import numpy as np
import subprocess


in_dir = "/data5/galaxy/project/mettl3_enrich/result_parse/bedtools_fisher"
sig_list, con_list, ratio_list, name_list = [], [], [], []
for tf_file in glob.glob("%s/*.bed" % in_dir):
    name = os.path.basename(tf_file).split(".bed")[0]
    command = "tail -1 %s" % tf_file
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    last_line = sub_p.communicate()[0].split()
    if float(last_line[2]) < 0.05:
        sig_list.append(os.path.basename(tf_file).split(".")[0])
        ratio_list.append(float(last_line[3]))
        name_list.append(name)
    else:
        con_list.append(os.path.basename(tf_file).split(".")[0])
print("significant: %d\ncontrol: %d" % (len(sig_list), len(con_list)))
index_list = np.argsort(ratio_list)
largest_index_list = index_list[-20:]
for i_index in largest_index_list:
    print("%s\t%f" % (name_list[i_index], ratio_list[i_index]))
