#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import subprocess
# from scipy import stats
from multiprocessing import Pool

##############################################################################################
#############################
tfbs_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_hg19_to_GRCh38"
mettle3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
tfbs_list = glob.glob("%s/*.bed" % tfbs_dir)
result_file = "/data5/galaxy/project/mettl3_enrich/result_parse/statistic_result.txt"
##############################################################################################


def statistic_line_num(infile):
    with open(infile, 'r') as f:
        line_num = len(f.readlines())
    return line_num


def statistic_variables(tfbs_bed):
    total_tfbs_num = statistic_line_num(tfbs_bed)
    num = subprocess.check_output("bedtools intersect -a %s -b %s -wa | wc -l" % (tfbs_bed, mettle3_bed), shell=True)
    proportion_tf = float(num) / total_tfbs_num
    return "%s\t%f\n" % (os.path.basename(tfbs_bed).split(".bed")[0], proportion_tf)


pool, result_list = Pool(), []
for tfbs in tfbs_list:
    result = pool.apply_async(statistic_variables, (tfbs, ))
    # ratio = statistic_variables(tfbs)
    result_list.append(result)
    # print("%s\t%f\n" % (os.path.basename(tfbs).split(".bed")[0], ratio))
pool.close()
pool.join()

with open(result_file, 'w') as fw:
    results = [x.get() for x in result_list]
    sorted_results = sorted(results, key=lambda d: float(d.split("\t")[1]))
    for x in sorted_results:
        fw.write(x.strip() + "\n")
