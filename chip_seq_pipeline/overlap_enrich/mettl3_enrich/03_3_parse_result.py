#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import subprocess
import numpy as np

in_dir = "/data5/galaxy/project/mettl3_enrich/result_parse/bedtools_jaccard"
complex_dir = "/data5/galaxy/project/what/TF_narrowPeak/clean_data"
complex_gene_list = [os.path.basename(x).split(".bed")[0] for x in glob.glob("%s/*.bed" % complex_dir)]
print(complex_gene_list[0])

complex_list, control_list = [], []
for tf_file in glob.glob("%s/*.bed" % in_dir):
    name = os.path.basename(tf_file).split(".bed")[0].upper()
    print(name)
    command = "tail -1 %s" % tf_file
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    ratio = float(sub_p.communicate()[0].split()[2])
    if name in complex_gene_list:
        complex_list.append(ratio)
    else:
        control_list.append(ratio)
print(np.average(complex_list), np.average(control_list))
