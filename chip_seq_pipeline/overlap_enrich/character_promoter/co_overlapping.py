#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import subprocess
from scipy import stats


reference_genome = "/data/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
total_promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/promoter_bed"
m6a_promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_gene/promoter_bed"
promoter_mettl3_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/promoter_mettl3"
tissue_list = ["brain", "heart", "kidney", "liver", "lung", "muscle", "placenta", "stomach"]
if not os.path.exists(promoter_mettl3_dir):
    os.makedirs(promoter_mettl3_dir)


def identify_co_overlapping():
    for tissue in tissue_list:
        promoter_mettl3_bed = get_occupied_promoter(tissue)
        str_promoter_m6a = get_m6a_promoter(tissue)
        str_promoter_free = get_free_promoter(tissue)
        a, c = stat_overlap_num(str_promoter_m6a, promoter_mettl3_bed)
        b, d = stat_overlap_num(str_promoter_free, promoter_mettl3_bed)
        print(a, b, c, d)
        oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
        print("%s\t%f\t%f\n" % (tissue, oddsratio, pvalue))


def get_occupied_promoter(tissue_name):
    total_promoter_bed = os.path.join(total_promoter_dir, "%s.bed" % tissue_name)
    command_pos = "bedtools intersect -a %s -b %s -wa | sort | uniq " % (total_promoter_bed, mettl3_bed)
    sub_p = subprocess.Popen(command_pos, shell=True, stdout=subprocess.PIPE)
    str_promoter_mettl3 = str(sub_p.communicate()[0])
    promoter_mettl3_bed = os.path.join(promoter_mettl3_dir, os.path.basename(total_promoter_bed))
    # print("########################################")
    # print(str_promoter_mettl3)
    # print("########################################")
    with open(promoter_mettl3_bed, 'w') as fw:
        fw.write(str_promoter_mettl3)
    return promoter_mettl3_bed


def get_m6a_promoter(tissue_name):
    m6a_promoter_bed = os.path.join(m6a_promoter_dir, "%s.bed" % tissue_name)
    print(m6a_promoter_bed)
    with open(m6a_promoter_bed, 'r') as f:
        print(m6a_promoter_bed)
        str_promoter_m6a = f.read()
        # print(str_promoter_m6a)
        # print("########################################")
        print("a+c = %d" % (int(len(str_promoter_m6a.strip().split("\n")))))
    return str_promoter_m6a


def get_free_promoter(tissue_name):
    total_promoter_bed = os.path.join(total_promoter_dir, "%s.bed" % tissue_name)
    # print(total_promoter_bed)
    m6a_promoter_bed = os.path.join(m6a_promoter_dir, "%s.bed" % tissue_name)
    # command = "bedtools intersect -v -a %s -b %s | sort | uniq " % (total_promoter_bed, m6a_promoter_bed)
    # sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    # str_promoter_free = sub_p.communicate()[0]
    with open(total_promoter_bed, 'r') as f_t:
        total_promoters = f_t.readlines()
        with open(m6a_promoter_bed, 'r') as f_m:
            m6a_promoters = f_m.readlines()
    # free_promoters = [x.strip() for x in total_promoters if x not in m6a_promoters]
    free_promoters = list(set(total_promoters).difference(set(m6a_promoters)))
    str_promoter_free = "\n".join([x.strip() for x in free_promoters])
    # print("########################################")
    # print(str_promoter_free)
    # print("########################################")
    print("b+d = %d" % (int(len(str_promoter_free.strip().split("\n")))))
    return str_promoter_free


def stat_overlap_num(str_total_bed, promoter_mettl3_bed):
    command = "bedtools intersect -a stdin -b %s -wa | sort | uniq | wc -l" % promoter_mettl3_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    overlap_num = int(sub_p.communicate(str_total_bed)[0].strip())
    return overlap_num, (len(str_total_bed.strip().split("\n")) - overlap_num)


if __name__ == '__main__':
    identify_co_overlapping()
