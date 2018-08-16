#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

# http://www.bio-info-trainee.com/1211.html
import os
import glob


peak = "/data5/galaxy/project/GWAS_analysis/data/ip/autosomal_peak/all_tissue/all_tissues.bed"
gwas_snp_bed="/data4/galaxy/project/new_analysis_m6a/gwas_analysis/snp_classification/merge/gwas_snp_total_withName.bed"
result_dir = "/data5/galaxy/project/tf_pathway/m6a_snp_intersect"
result_bed_dir = "%s/peak" % result_dir
result_snp_dir = "%s/snp" % result_dir
for i_dir in [result_dir, result_bed_dir, result_snp_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)


# remove file duplicate rows
def intersect_peak(peak_bed):
    result_file = os.path.join(result_bed_dir, os.path.basename(peak_bed))
    os.system("bedtools intersect -a %s -b %s -wa | awk '!seen[$0]++' - > %s" % (peak_bed, gwas_snp_bed, result_file))


def intersect_snp(peak_bed):
    result_file = os.path.join(result_snp_dir, os.path.basename(peak_bed))
    os.system("bedtools intersect -a %s -b %s -wa | awk '!seen[$0]++' - > %s" % (gwas_snp_bed, peak_bed, result_file))


intersect_peak(peak)
intersect_snp(peak)
#