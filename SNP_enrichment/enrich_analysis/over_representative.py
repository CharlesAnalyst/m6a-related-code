#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import pandas as pd
from scipy.stats import binom_test


# #######function########
# test over representation of SNPs with respect to genomic background.
# determine the genomic space spanned by the m6a regions.
# use the control regions as genomic background.
#######################
peak_bed = "/data5/galaxy/project/GWAS_unfilter/data/ip/autosomal_peak/all_tissue/all_tissues.bed"
control_bed = "/data5/galaxy/project/bam2bed_input/merged_data/all_tissue/all_tissues.bed"
autosomal_snp_bed = "/data/database/snp/autosomal_GRCh38_snp.bed"
result_dir = "/data5/galaxy/project/snp_analysis/binomial_test"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = os.path.join(result_dir, "binomial_test_result_2.txt")
##################################################################
# cat consecutive_bed/*.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > all_tissue/all_tissues.bed


def stat_bed_length(i_bed):
    df = pd.read_table(i_bed, sep="\t", header=None, names=["chr", "start", "end"])
    df["length"] = df["end"] - df["start"]
    total_length = df["length"].sum()
    return total_length


def stat_overlap_num():
    control_popen = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (autosomal_snp_bed, control_bed))
    peak_popen = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (autosomal_snp_bed, peak_bed))
    control_num, peak_num = int(control_popen.read().split()[0]), int(peak_popen.read().split()[0])
    return control_num, peak_num


def binomial_test(success_num, total_num, prob):
    p_value = binom_test(success_num, total_num, prob, alternative='greater')
    return p_value


def main_run():
    peak_len, control_len = stat_bed_length(peak_bed), stat_bed_length(control_bed)
    len_ratio = float(peak_len) / float(control_len)
    control_num, peak_num = stat_overlap_num()
    p_value = binomial_test(peak_num, control_num, len_ratio)
    result_line = "%d\t%d\t%f\t%f" % (peak_num, control_num, len_ratio, p_value)
    return result_line


if __name__ == "__main__":
    result = main_run()
    with open(result_file, 'w') as fw:
        fw.write("binomial test:\n")
        fw.write("peak_num\tcontrol_num\tlen_ratio\tp_value\n")
        fw.write(result)


