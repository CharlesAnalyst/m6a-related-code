#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os

#work_dir = "/data5/galaxy/project/GWAS_analysis/universal/Fisher_exact_test"
work_dir = "/data5/galaxy/project/GWAS_analysis/cluster_by_map/Fisher_exact_test"
os.chdir(work_dir)
infile = "Fisher_m6a_vs_input.txt"
outfile = "Fisher_m6a_vs_input_matrix.txt"

with open(infile, 'r') as f:
    f.readline()
    contents = f.readlines()
    contents.sort(key=lambda d: (d.split("\t")[0], d.split("\t")[1]))
    tissue_list, trait_list, oddsratio_list = [], [], []
    for line in contents:
        info = line.strip().split("\t")
        tissue_list.append(info[0])
        trait_list.append(info[1])
        oddsratio_list.append(info[2])
    tissue_list = list(set(tissue_list))
    print("oddsratio_list length: %d" % len(oddsratio_list))
    tissue_list.sort()
    trait_list = list(set(trait_list))
    trait_list.sort()

    with open(outfile, 'w') as fw:
        fw.write("tissue\t%s\n" % ("\t".join(trait_list)))
        for i in range(len(tissue_list)):
            i_oddsratio_list = [oddsratio_list[j] for j in range(i * len(trait_list), (i+1)*len(trait_list))]
            oddsratio_value = "\t".join(i_oddsratio_list)
            fw.write("%s\t%s\n" % (tissue_list[i], oddsratio_value))


