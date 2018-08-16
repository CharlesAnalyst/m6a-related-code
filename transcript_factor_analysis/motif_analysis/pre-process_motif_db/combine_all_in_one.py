#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
os.chdir("D:\\Project\\Motif\\new_download_data\\JASPAR\\combine_result_motifs_clean")
out_file = "D:\\Project\\Motif\\new_download_data\\JASPAR\\all_in_one\\motifs.txt"


in_files = glob.glob("*.txt")
with open(out_file, 'a') as fw:
    with open(in_files[0], 'r') as f:
        [fw.write(line) for line in f.readlines()]
    for i in range(1, len(in_files)):
        with open(in_files[i], 'r') as f:
            contents = f.readlines()
            for j in range(9, len(contents)):
                fw.write(contents[j] + "\r")
            fw.write("\r")

