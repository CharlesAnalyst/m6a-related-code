#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import datetime
import os


filtered_file = "/data5/galaxy/project/GWAS_db/replicated_GWS_snp.txt"
result_dir = "/data5/galaxy/project/GWAS_db/snp_classification_2"


def read_disease_file():
    disease_list = []
    with open(filtered_file, 'r') as f:
        f.readline()
        for line in f.readlines():
            info = line.strip().split("\t")
            chr_id, chr_pos, rs_id, mapped_trait = ("chr%s" % info[11]), info[12], info[21], info[34].lower()
            chr_start, chr_end = chr_pos, str(int(chr_pos) + 1)     # GWAS Catalog us
            disease_list.append("%s\t%s\t%s\t%s\t%s" % (chr_id, chr_start, chr_end, rs_id, mapped_trait))
    return disease_list