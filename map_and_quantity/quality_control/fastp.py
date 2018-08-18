#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool


##################################################################################
fq_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/m6a_level/fq"     # "SRR1596087.fastq.gz"
postfix = ".fastq.gz"
fq_list = glob.glob("%s/*%s" % (fq_dir, postfix))
result_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/m6a_level/clean_fq"
thread_num = 40 / len(fq_list)
#######################################


def get_prefix(fq):
    prefix = os.path.basename(fq).split(postfix)[0]
    return prefix


def run_fastp(fq):
    prefix = get_prefix(fq)
    sub_dir = os.path.join(result_dir, prefix)
    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)
    out_name = os.path.join(sub_dir, prefix)
    filtered_fq, html_report, json_report = "%s.fastq" % out_name, "%s.html" % out_name, "%s.json" % out_name
    os.system("fastp -i %s -o %s -h %s -j %s -w %d -z 2" % (fq, filtered_fq, html_report, json_report, thread_num))


if __name__ == "__main__":
    pool = Pool(processes=len(fq_list))
    for i_fq in fq_list:
        pool.apply_async(run_fastp, (i_fq, ))
    pool.close()
    pool.join()
    # run_fastp("/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/m6a_level/fq/SRR1596097.fastq.gz")


