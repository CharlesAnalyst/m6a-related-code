#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob

# total
original_gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/gene_bed"
original_promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/promoter_bed"
# m6a
gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_gene/gene_bed"
promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_gene/promoter_bed"
# free
free_gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_free/gene_bed"
free_promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_free/promoter_bed"
for i_dir in [free_gene_dir, free_promoter_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)


gene_dir_list = [original_gene_dir, gene_dir]
promoter_dir_list = [original_promoter_dir, promoter_dir]


def generate_free_data():
    gene_dict_list = [generate_dict(g_dir) for g_dir in gene_dir_list]
    promoter_dict_list = [generate_dict(pro_dir) for pro_dir in promoter_dir_list]
    a = [gene_dict_list, promoter_dict_list]
    b = [free_gene_dir, free_promoter_dir]
    for i in range(len(a)):
        process_one_type_data(a[i], b[i])


def process_one_type_data(data_dict_list, result_dir):
    total_dict, m6a_dict = data_dict_list[0], data_dict_list[1]
    for tissue in total_dict:
        total_bed, m6a_bed = total_dict[tissue], m6a_dict[tissue]
        free_bed = os.path.join(result_dir, os.path.basename(total_bed))
        get_difference_set(total_bed, m6a_bed, free_bed)


def get_difference_set(total_bed, m6a_bed, free_bed):
    with open(total_bed, 'r') as f_total:
        total_list = f_total.readlines()
        with open(m6a_bed, 'r') as f_m6a:
            m6a_list = f_m6a.readlines()
            free_list = list(set(total_list).difference(set(m6a_list)))
    free_list.sort(key=lambda d: (d.split("\t")[0].split("chr")[1], d.split("\t")[1]))
    with open(free_bed, 'w') as fw:
        fw.writelines(free_list)


def generate_dict(data_dir):
    result_dict = {}
    gene_list = glob.glob("%s/*.bed" % data_dir)
    for g_bed in gene_list:
        result_dict[os.path.basename(g_bed).split(".bed")[0].lower()] = g_bed
    return result_dict


if __name__ == '__main__':
    generate_free_data()