#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import subprocess
from multiprocessing import Pool


###################################################################################################
data_type = "total_gene"  # m6a_gene
#######################
gene_dir, promoter_dir, result_dir = "", "", ""
tfbs_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean"
tfbs_list = glob.glob("%s/*.bed" % tfbs_dir)
if data_type == "m6a_gene":
    gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/random_m6a/gene_bed"
    promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/random_m6a/promoter_bed"
    result_dir = "/data5/galaxy/project/promoter_TF_enrich/comprehensive_enrich/%s" % data_type   # ####
elif data_type == "total_gene":
    gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/random_total/gene_bed"
    promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/random_total/promoter_bed"
    result_dir = "/data5/galaxy/project/promoter_TF_enrich/comprehensive_enrich/%s" % data_type
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
######################################################################################################


def match_file_according_tissue(tissue_name, promoter_dict):
    promoter_bed = ""
    try:
        promoter_bed = promoter_dict[tissue_name]
    except KeyError:
        print("%s not in promoter_dict!" % tissue_name)
    return promoter_bed


def get_dict(data_dir):
    tissue_dict, file_list = {}, glob.glob("%s/*.bed" % data_dir)
    for i_file in file_list:
        tissue_name = os.path.basename(i_file).split(".bed")[0].lower()
        tissue_dict[tissue_name] = i_file
    return tissue_dict


def stat_overlap_number(promoter_bed):
    tissue_name, result_list = os.path.basename(promoter_bed).split(".bed")[0].lower(), []
    result_file = os.path.join(result_dir, "overlap-number_%s.txt" % tissue_name)
    if not os.path.exists(result_file):
        with open(promoter_bed, 'r') as f:
            for line in f.readlines():
                result_list.append(stat_each_gene(line.strip()))
        write_to_file(result_list, result_file)


def stat_each_gene(each_gene_line):
    gene, overlap_num, pool, result_list = each_gene_line.split("\t")[3], 0, Pool(), []
    for tfbs_bed in tfbs_list:
        score = pool.apply_async(stat_each_tfbs, (each_gene_line, tfbs_bed))
        result_list.append(score)
    pool.close()
    pool.join()
    sum_score = sum([int(x.get()) for x in result_list])
    return "%s\t%d\n" % (gene, sum_score)


def stat_each_tfbs(each_gene_line, tfbs_bed):
    score = 0
    command = "bedtools intersect -a %s -b stdin -f 0.5 -wa | sort | uniq | wc -l" % tfbs_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    tf_num = int(sub_p.communicate(each_gene_line)[0])
    if tf_num > 0:
        score = 1
    return score


def write_to_file(result_list, result_file):
    with open(result_file, 'w') as fw:
        fw.writelines(result_list)


if __name__ == '__main__':
    tissue_list = [os.path.basename(x).split(".bed")[0].lower() for x in glob.glob("%s/*.bed" % gene_dir)]
    pro_dict = get_dict(promoter_dir)
    for tissue in tissue_list:
        print(tissue)
        promoter_file = match_file_according_tissue(tissue, pro_dict)
        stat_overlap_number(promoter_file)
