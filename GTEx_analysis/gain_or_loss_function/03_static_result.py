#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import re
import glob


data_file = "/data5/galaxy/project/snp_analysis/GTEx_analysis/gain_or_loss/02_parse_result/GRCh38.bed"

group_dict = {}
with open(data_file, 'r') as f:
    for line in f.readlines():
        info_list = re.split(r'[;\t]', line.strip())
        tissue = info_list[3]
        function_t = info_list[-2]
        i_key = "%s-%s" % (tissue, function_t)
        group_dict[i_key] = group_dict.get(i_key, 0) + 1

# for a, b in group_dict.items():
    # print(a, b)
# print("################################")

result_list = []
for i_key, count in group_dict.items():
    result_list.append("%s\t%d" % (i_key, count))

result_list.sort()
result_dict = {}
term_list = list(set([x.split("-Gain")[0].split("-Loss")[0] for x in result_list]))
for term in term_list:
    result_dict[term] = {}
    if "%s-Gain" % term in group_dict:
        result_dict[term]["gain"] = group_dict["%s-Gain" % term]
    else:
        result_dict[term]["gain"] = 0
    if "%s-Loss" % term in group_dict:
        result_dict[term]["loss"] = group_dict["%s-Loss" % term]
    else:
        result_dict[term]["loss"] = 0

# print(result_dict)
for i_key, i_dict in result_dict.items():
    print("%s\t%d\t%d" % (i_key, i_dict["gain"], i_dict["loss"]))
