#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
from scipy import stats
from multiprocessing import Pool


# #######function########
# test whether GWAS SNPs are overrepresented within dynamic regions compared to HapMap SNPs: binomial test b(x;n,p)
# expectation p equal to the fraction of autosomal HapMap SNPs present within our full dynamics region set;
#########################
peak_bed = "/data5/galaxy/project/data/total_m6a_peak/universal/universal.bed"
gwas_snp_bed = "/data5/galaxy/project/GWAS_db/03_replicated_GRCh38/autosomol_gwas_snp.bed"
autosomal_snp_bed = "/data/database/snp/autosomal_GRCh38_snp.bed"
total_autosomal_num, total_gwas_num = 81039058, 9760
# result_dir = "/data5/galaxy/project/GWAS_analysis/fisher_exact_test/gwas_vs_snp"
# if not os.path.exists(result_dir):
#     os.makedirs(result_dir)
# result_file = os.path.join(result_dir, "Fisher_exact_test_result.txt")
##################################################################
# cat gwas_snp_*.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > merge/gwas_snp_total.bed


def stat_overlap_num():
    gwas_popen = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (gwas_snp_bed, peak_bed))
    autosomal_popen = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (autosomal_snp_bed, peak_bed))
    gwas_num, auto_num = int(gwas_popen.read().split()[0]), int(autosomal_popen.read().split()[0])
    print(gwas_num, auto_num)
    return gwas_num, auto_num


def fisher_exact_test(gwas_num, auto_num):
    #     gwas_snp   snp
    # m6a    a        b
    # con    c        d
    print([[gwas_num, auto_num], [total_gwas_num - gwas_num, total_autosomal_num - auto_num]])
    oddsratio, pvalue = stats.fisher_exact([[gwas_num, auto_num], [total_gwas_num - gwas_num, total_autosomal_num - auto_num]])
    # return result_line
    print(oddsratio, pvalue)


def main_method():
    gwas_num, auto_num = stat_overlap_num()
    fisher_exact_test(gwas_num, auto_num)


if __name__ == "__main__":
    main_method()
