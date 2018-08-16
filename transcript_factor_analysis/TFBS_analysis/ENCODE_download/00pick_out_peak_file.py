#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os

#########################################################
metadata_csv = "E:\\Users\\cshuo\\Desktop\\metadata.tsv"
result_file = "E:\\Users\\cshuo\\Desktop\\get_target.txt"
#########################################################
relation_dict, target_expre_dict, expre_rep_url_dict, expre_url_dict = {}, {}, {}, {}
with open(metadata_csv, 'r') as f:
    title = f.readline()
    contents = f.readlines()
    for line in contents:
        info = line.split("\t")
        expre_target, expre_acc, out_type, assembly, bio_rep, url = info[4], info[2], info[1], info[5], info[6], info[8]
        target_expre_dict[expre_target] = target_expre_dict.get(expre_target, []) + [expre_acc]
        if out_type == "peaks" and assembly == "hg19":
            expre_rep_url_dict[expre_acc] = expre_rep_url_dict.get(expre_acc, []) + ["%s\t%s" % (bio_rep, url)]
# remove duplicate
for i_key, i_values in target_expre_dict.items():
    target_expre_dict[i_key] = list(set(i_values))
# get url
# assert "TAF1-human" in target_expre_dict
# assert "ENCSR000AKO" in expre_rep_url_dict
for expre, bio_rep_url_list in expre_rep_url_dict.items():
    if len(bio_rep_url_list) <= 1:
        expre_url_dict[expre] = bio_rep_url_list[0].split("\t")[1]
    else:
        for i in range(len(bio_rep_url_list)):
            if len(bio_rep_url_list[i].split("\t")[0].split(",")) >= 2:
                expre_url_dict[expre] = bio_rep_url_list[i].split("\t")[1]
#
with open(result_file, 'w') as fw:
    for expre_target, expre_list in target_expre_dict.items():
        for expre in expre_list:
            if expre in expre_url_dict:
                fw.write("%s\t%s\t%s\n" % (expre_target, expre, expre_url_dict[expre]))
            else:
                # These experiment accession can't pick out peak bed based on above method!
                print(expre)
