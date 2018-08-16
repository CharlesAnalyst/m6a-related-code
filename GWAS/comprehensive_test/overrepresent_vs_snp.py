#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
from multiprocessing import Pool
from scipy.stats import binom_test


# #######function########
# test whether GWAS SNPs are overrepresented within dynamic regions compared to HapMap SNPs: binomial test b(x;n,p)
# expectation p equal to the fraction of autosomal HapMap SNPs present within our full dynamics region set;
# n equal to the number of autosomal GWAS SNPs;
# x equal to GWAS SNPs present in dynamic regions.
#########################
peak_dir = "/data5/galaxy/project/data/total_m6a_peak/autosomal_peak"
gwas_snp_bed = "/data5/galaxy/project/GWAS_db/03_replicated_GRCh38/all_replicated_gwas_snp.bed"
autosomal_snp_bed = "/data/database/snp/autosomal_GRCh38_snp.bed"
total_autosomal_num, total_gwas_num = 81039058, 9823
result_dir = "/data5/galaxy/project/GWAS_analysis/binomial_test/gwas_vs_snp"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = os.path.join(result_dir, "binomial_test_result_2.txt")
##################################################################


def stat_overlap_num(peak_bed):
    gwas_popen = os.popen("bedtools intersect -a %s -b %s -wa |sort|uniq| wc -l" % (gwas_snp_bed, peak_bed))
    autosomal_popen = os.popen("bedtools intersect -a %s -b %s -wa |sort|uniq| wc -l" % (autosomal_snp_bed, peak_bed))
    gwas_num, auto_num = int(gwas_popen.read().split()[0]), int(autosomal_popen.read().split()[0])
    return gwas_num, auto_num


def binomial_test(gwas_num, auto_num):
    p = float(auto_num) / float(total_autosomal_num)
    print(auto_num, total_autosomal_num, p)
    p_value = binom_test(gwas_num, total_gwas_num, p, alternative='greater')
    result_line = "%d\t%d\t%f\t%f" % (gwas_num, total_gwas_num, p, p_value)
    return result_line


def process_each_tissue(peak_bed):
    gwas_num, auto_num = stat_overlap_num(peak_bed)
    tissue_name = os.path.basename(peak_bed).split(".bed")[0]
    result_line = binomial_test(gwas_num, auto_num)
    return "%s\t%s\n" % (tissue_name, result_line)


if __name__ == "__main__":
    pool, results, name_list = Pool(processes=20), [], []
    bed_list = glob.glob("%s/*.bed" % peak_dir)
    results = pool.map_async(process_each_tissue, bed_list)
    pool.close()
    pool.join()
    result_list = results.get()
    with open(result_file, 'w') as fw:
        fw.write("tissue\tgwas_num\ttotal_gwas_num\tprobability(autosomal)\tp_value\n")
        fw.writelines(result_list)
