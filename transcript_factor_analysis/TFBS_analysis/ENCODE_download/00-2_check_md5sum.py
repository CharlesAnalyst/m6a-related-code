#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import pandas as pd


in_file = "/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/get_target.txt"
metadata_csv = "/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/metadata.tsv"
md5sum_file = "/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/raw_data/md5_sum.txt"


df = pd.read_table(in_file, sep="\t", header=None, names=["1", "2", "url"])
url_list = df.url.tolist()

url_md5_dict = {}
with open(metadata_csv, 'r') as f:
    title = f.readline()
    contents = f.readlines()
    for line in contents:
        info = line.split("\t")
        url, md5 = info[8].strip(), info[37].strip()
        url_md5_dict[url] = md5

#
with open(md5sum_file, 'w') as fw:
    for url in url_list:
        md5 = url_md5_dict[url]
        file = url.split("/")[-1].strip().strip("$")
        print(file)
        fw.write("%s %s\n" % (md5, file))

