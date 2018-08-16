#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import re
import math
import subprocess
import numpy as np
from scipy.stats import ttest_ind


reference_genome = "/data/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
promoter_bed = "/data5/galaxy/project/data/promoter/2k_100/genes_promoters.bed"
mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
result_dir = "/data5/galaxy/project/promoter_TF_enrich/mettl3_vs_free"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = os.path.join(result_dir, "mettl3_vs_free.txt")


def compare_CpG_enrichment_in_promoter():
    score_list = []
    str_promoter_pos, str_promoter_neg = get_occupied_promoter()
    for str_promoter in [str_promoter_pos, str_promoter_neg]:
        score_list.append(enrich_CpG_in_promoter(str_promoter))
    write_data_into_file(score_list[0], score_list[1])
    t_statistic, p_value = ttest_ind(score_list[0], score_list[1])
    print("#############################################\n")
    print("promoter occupied by METTL3:\t%f" % np.mean(score_list[0]))
    print("promoter not occupied by METTL3:\t%f" % np.mean(score_list[1]))
    print("t test:")
    print(t_statistic, p_value)


def enrich_CpG_in_promoter(str_promoter):
    str_seq_list, score_list = get_sequence_from_bed(str_promoter).split("\n"), []
    for i in range(0, len(str_seq_list)-1, 2):
        title, seq = str_seq_list[i], str_seq_list[i+1]
        score_list.append(calculate_CpG_density(seq))
    return score_list


def get_occupied_promoter():
    command_pos = "bedtools intersect -a %s -b %s -wa | sort | uniq " % (promoter_bed, mettl3_bed)
    command_neg = "bedtools intersect -a %s -b %s -v | sort | uniq " % (promoter_bed, mettl3_bed)
    sub_p = subprocess.Popen(command_pos, shell=True, stdout=subprocess.PIPE)
    str_promoter_pos = sub_p.communicate()[0]
    sub_p = subprocess.Popen(command_neg, shell=True, stdout=subprocess.PIPE)
    str_promoter_neg = sub_p.communicate()[0]
    return str_promoter_pos, str_promoter_neg


def get_sequence_from_bed(str_promoter):
    command = "bedtools getfasta -fi %s -bed stdin" % reference_genome
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    str_sequence = sub_p.communicate(str_promoter)[0]
    return str_sequence


def calculate_CpG_density(sequence):
    seq = sequence.lower()
    cg_num, c_num, g_num = len(re.findall("cg", seq)), len(re.findall("c", seq)), len(re.findall("g", seq))
    score = cg_num / (math.pow((c_num + g_num) / 2.0, 2) / (len(seq) * 1.0))
    return score


def write_data_into_file(pos_scores, neg_scores):
    with open(result_file, 'w') as fw:
        for score in pos_scores:
            fw.write("mettl3\t%f\n" % score)
        for score in neg_scores:
            fw.write("mettl3 free\t%f\n" % score)


if __name__ == '__main__':
    compare_CpG_enrichment_in_promoter()
