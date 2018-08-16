#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os

total_motif = "/data/database/tf_motif/motif_matrix.txt"
part_motif = "/data5/galaxy/project/tf_analysis/interacted_TF/part_motif.txt"

tf_list = ["THAP11", "NRF1", "ETV6", "RELA", "TP53", "MAFB", "OTX1", "STAT1", "STAT4", "MYC", "WT1", "JUN", "IKZF1"]
with open(part_motif, 'w') as fw:
    with open(total_motif, 'r') as f:
        contents = f.readlines()
        # write title line
        for i in range(9):
            fw.write(contents[i])
        j = 0
        while j < len(contents):
            if contents[j].startswith("MOTIF") and contents[j].split("\t")[1].strip() in tf_list:
                fw.write(contents[j])
                j += 1
                while not contents[j].startswith("MOTIF"):
                    fw.write(contents[j])
                    j += 1
            j += 1