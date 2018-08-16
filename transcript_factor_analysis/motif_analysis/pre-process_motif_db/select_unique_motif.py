#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import shutil

motif_dir = "D:\\Project\\Motif\\CIS_BP_database\\pm\\pwms_all_motifs"
result_motif_dir = "D:\\Project\\Motif\\CIS_BP_database\\pm\\pwms_unique_motifs"
if not os.path.exists(result_motif_dir):
    os.makedirs(result_motif_dir)
tf_info = "D:\\Project\\Motif\\CIS_BP_database\\pm\\TF_Information.txt"
unique_motif = "D:\\Project\\Motif\\CIS_BP_database\\pm\\unique_TF_Information.txt"

relation_dict = {}
with open(tf_info, 'r') as f:
    first_list = []
    for line in f.readlines():
        info = line.strip().split("\t")
        tf_id, motif_id, tf_name, tf_status, motif_type, source_id = info[0], info[3], info[6], info[8], info[14], info[15]
        if source_id in ["Transfac", "HocoMoco", "JASPAR", "ENCODE", "modENCODE"]:
            if tf_status == "D":
                new_line = "\t".join([tf_id, motif_id, tf_name, motif_type, source_id])
                relation_dict[tf_name] = relation_dict.get(tf_name, []) + [new_line]


result_list = []
for tf_name, lines in relation_dict.items():
    if len(lines) == 1:
        result_list.append(lines[0])
    else:
        transfac_list, jaspar_list, encode_list, modencode_list, hocomoco_list = [], [], [], [], []
        for line in lines:
            tf_id, motif_id, tf_name, motif_type, source_id = line.split("\t")
            if motif_type == "Transfac":
                transfac_list.append(line)
            if motif_type == "HocoMoco":
                hocomoco_list.append(line)
            if motif_type == "JASPAR":
                jaspar_list.append(line)
            if motif_type == "ENCODE":
                encode_list.append(line)
            if motif_type == "modENCODE":
                modencode_list.append(line)
        if len(transfac_list) > 0:
            result_list.append(transfac_list[0])
        elif len(jaspar_list) > 0:
            result_list.append(jaspar_list[0])
        elif len(encode_list) > 0:
            result_list.append(encode_list[0])
        elif len(modencode_list) > 0:
            result_list.append(modencode_list[0])
        elif len(hocomoco_list) > 0:
            result_list.append(hocomoco_list[0])


with open(unique_motif, 'w') as fw:
    fw.write("TF_ID\tMotif_ID\tTF_Name\tMotif_Type\tMSource_Identifier\n")
    for line in result_list:
        fw.write(line + "\n")

os.chdir(motif_dir)
# motif_files = [os.path.abspath(file) for file in glob.glob(os.path.join(motif_dir, "M*.txt"))]
motif_files = glob.glob("M*.txt")
for line in result_list:
    motif_file_name = "%s.txt" % line.split("\t")[1]
    if motif_file_name in motif_files:
        shutil.copyfile(motif_file_name, os.path.join(result_motif_dir, motif_file_name))
    else:
        print("%s file is not in motif directory!" % motif_file_name)
