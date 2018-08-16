#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import shutil

jaspar_dir = "D:\\Project\\Motif\\new_download_data\\JASPAR\\combine_result_motifs"
os.chdir(jaspar_dir)
result_dir = "D:\\Project\\Motif\\new_download_data\\JASPAR\\combine_result_motifs_clean"
header_string = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n\nA 0.25 C 0.25 G 0.25 T " \
                "0.25\n\n"

files = glob.glob("*.txt")
count = 1
for file in files:
    motif_dict = {}
    with open(file, 'r') as f:
        contents = f.readlines()
        for i in range(len(contents)):
            if contents[i].startswith("MOTIF"):
                motif_name = contents[i].split(" ")[1].strip()
                motif_dict[motif_name] = motif_dict.get(motif_name, []) + [contents[i].strip()]
                i += 1
                while i < len(contents) and not contents[i].startswith("MOTIF"):
                    motif_dict[motif_name] += [contents[i].strip()]
                    i += 1
    unique_key, unique_value = "", []
    for i_key, i_value in motif_dict.items():
        if "full" in i_key:
            unique_key, unique_value = i_key, i_value
    if unique_key == "":
        for i_key, i_value in motif_dict.items():
            if "DBD" in i_key:
                unique_key, unique_value = i_key, i_value
    # force class
    if unique_key == "":
        for i_key, i_value in motif_dict.items():
            if ".1" in i_key:
                unique_key, unique_value = i_key, i_value
    # classed
    if unique_key != "":
        count += 1
        out_file = os.path.join(result_dir, file)
        with open(out_file, 'w') as fw:
            fw.write(header_string + "\n")
            for line in unique_value:
                fw.write(line + "\n")
    elif len(list(motif_dict.keys())) == 1:
        count += 1
        out_file = os.path.join(result_dir, file)
        shutil.copyfile(file, out_file)
    else:
        out_file = os.path.join(result_dir, "unclass-%s" % file)
        shutil.copyfile(file, out_file)
print(count)
