#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import itertools

relation_file = "relation_repeat.txt"
bed_dir = "/data5/galaxy/project/TF_narrowPeak/repeat"
os.chdir(bed_dir)
result_dir = bed_dir

rel_dict = {}
with open(relation_file, 'r') as f:
    for line in f.readlines():
        info = line.strip().split("\t")
        rel_dict[info[0]] = rel_dict.get(info[0], []) + ["%s.bed" % info[1]]


def multi_intersect(tissue_name, bed_list):
    total_list = list(itertools.permutations(bed_list, 2))
    for i in range(len(total_list)):
        result_i = os.path.join(result_dir, "%s_%d.bed" % (tissue_name, i))
        os.system("intersectBed -a %s -b %s -f 0.5 > %s" % (total_list[i][0], total_list[i][1], result_i))
    result_file = os.path.join(result_dir, "%s.bed" % tissue_name)
    os.system("cat %s_*.bed | sort -k1,1 -k2,2n | mergeBed -i - > %s" % (tissue_name, result_file))
    os.system("rm %s_*.bed" % tissue_name)


for name, i_list in rel_dict.items():
    if len(i_list) == 1:
        result_file = "%s.bed" % name
        os.system("sort -k1,1 -k2,2n %s | mergeBed -i - > %s" % (i_list[0], result_file))
    else:
        multi_intersect(name, i_list)
