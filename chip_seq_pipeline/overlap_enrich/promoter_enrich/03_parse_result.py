#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd


work_dir = "/data5/galaxy/project/promoter_TF_enrich/tissue_enr_result"
file_list = glob.glob("%s/*.txt" % work_dir)

complex_dir = "/data5/galaxy/project/what/TF_narrowPeak/clean_data"
complex_gene_list = [os.path.basename(x).split(".bed")[0] for x in glob.glob("%s/*.bed" % complex_dir)]
print(complex_gene_list)

result_dir = "/data5/galaxy/project/promoter_TF_enrich/tissue_enr_result/two_class"
complex_dir, control_dir = "%s/complex" % result_dir, "%s/control" % result_dir
for i_dir in [complex_dir, control_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)

print("number of significant && non")
for x in file_list:
    complex_file, control_file = os.path.join(complex_dir, os.path.basename(x)), os.path.join(control_dir, os.path.basename(x))
    complex_list, control_list = [], []
    df = pd.read_table(x, sep="\t", header=None, names=["a", "b", "c", "d"])
    # df_sig, df_no = df[df["d"] <= 0.05], df[df["d"] > 0.05]
    # print("%s\t%d\t%d" % (os.path.basename(x).split(".")[0], len(df_sig), len(df_no)))
    df_a = df.set_index(["b"])
    with open(complex_file, 'w') as f_c:
        with open(control_file, 'w') as f_con:
            for name, values in df_a.iterrows():
                print(name)
                if name in complex_gene_list:
                    f_c.write("%s\t%s\n" % (name, "\t".join([str(m) for m in values])))
                else:
                    f_con.write("%s\t%s\n" % (name, "\t".join([str(m) for m in values])))
