#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd


def main_method():
    work_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/04_intersect"
    os.chdir(work_dir)
    tissue_list = ["brain", "heart", "kidney", "liver", "lung", "muscle", "placenta", "stomach"]
    files = []
    for tissue in tissue_list:
        files.append("%s/%s.txt" % (tissue, tissue))
    # files = glob.glob("%s/*_fimo.txt" % work_dir)
    # tissue_list = [os.path.basename(x).split("_")[0] for x in files]
    first_file = files[0]
    total_data = pd.read_table(first_file, sep="\t", header=None, names=["motif_id", "ratio"])
    title = ["motif_id", tissue_list[0]]
    for i in range(1, len(files)):
        title.append(tissue_list[i])
        i_data = pd.read_table(files[i], sep="\t", names=["motif_id", "ratio"])
        total_data = pd.merge(total_data, i_data, on="motif_id", how="outer")
    total_data = total_data.fillna(0.0)
    total_data.columns = title
    total_data.to_csv("bigTable.txt", sep="\t", index=False)


main_method()

