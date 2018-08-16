#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
# @function:
# 1. Transform different database id to protein name;
# 2.

import os
info_file = ""
##########
# os.chdir("D:\\Project\\Motif\\new_download_data\\HumanTF1")
# project = "HumanTF1"
# in_file = "total_motifs.txt"
##########
os.chdir("D:\\Project\\Motif\\new_download_data\\JASPAR\\new")
project = "JASPAR"
in_file = "all_JASPAR_MEME.txt"
##########
# os.chdir("D:\\Project\\Motif\\new_download_data\\uniprob")
# project = "UniPROBE"
# in_file = "total_motifs.txt"
##########
""""
os.chdir("D:\\Project\\Motif\\new_download_data\\MEME_HOCOMOCO")
project = "HOCOMOCO"
in_file = "HOCOMOCOv11_core_HUMAN_mono_meme_format.meme.txt"
info_file = "HUMAN_mono_motifs.tsv"
"""

out_dir = os.path.join("combine_result_motifs")
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
header_list, motif_dict = [], {}


with open(in_file, 'r') as f:
    """
    for i in range(10):
        header_list.append(f.readline().strip())
    """
    header_string = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n\nA 0.25 C " \
                    "0.25 G 0.25 T 0.25\n\n"
    header_list.append(header_string)
    contents = f.readlines()
    if project == "HumanTF1" or project == "JASPAR":
        for i in range(len(contents)):
            if contents[i].startswith("MOTIF"):
                protein_name = contents[i].split(" ")[-1].strip()
                motif_dict[protein_name] = motif_dict.get(protein_name, []) + [contents[i].strip()]
                i += 1
                while i < len(contents) and not contents[i].startswith("MOTIF"):
                    motif_dict[protein_name] += [contents[i].strip()]
                    i += 1
    elif project == "UniPROBE":
        for i in range(len(contents)):
            if contents[i].startswith("MOTIF") and "_REF_" in contents[i] and "_R1" in contents[i]:
                protein_name = contents[i].split(" ")[-1].strip()
                motif_dict[protein_name] = motif_dict.get(protein_name, []) + [contents[i].strip()]
                i += 1
                while i < len(contents) and not contents[i].startswith("MOTIF"):
                    motif_dict[protein_name] += [contents[i].strip()]
                    i += 1
    elif project == "HOCOMOCO":
        info_dict, contents_list = {}, []
        with open(info_file, 'r') as fr:
            for r_line in fr.readlines():
                model_name, tf_name = r_line.split("\t")[0], r_line.split("\t")[2]
                info_dict[model_name] = tf_name
        with open(in_file, 'r') as fr:
            contents_list = fr.readlines()
            for i in range(len(contents_list)):
                if contents_list[i].startswith("MOTIF"):
                    i_name = contents_list[i].split()[-1].strip()
                    if i_name in info_dict:
                        contents_list[i] = contents_list[i].replace(i_name, info_dict[i_name])
                    else:
                        print("%s doesn't have according information!" % contents_list[i])
        with open("motifs.txt", 'w') as f_w:
            for r_line in contents_list:
                f_w.write(r_line.strip() + "\n")


def write_header_and_content_to_file(out_file, out_list):
    with open(out_file, 'w') as fw:
        for line in header_list:
            fw.write(line)
        for line in out_list:
            fw.write(line.strip() + "\n")


os.chdir(out_dir)
for i_key, i_value in motif_dict.items():
    if "/" not in i_key:
        print(i_key)
        for protein in i_key.split("::"):
            out_put = "%s.txt" % protein
            write_header_and_content_to_file(out_put, i_value)
    else:
        print("un process protein name: %s" % i_key)



