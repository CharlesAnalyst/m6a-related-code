#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
from multiprocessing import Pool


tf_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean"
mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
chrom_size_file = "/data5/galaxy/project/mettl3_enrich/result_parse/hg38.chrom.sizes"
result_dir = "/data5/galaxy/project/mettl3_enrich/result_parse/bedtools_result"


def fisher_test(tf_bed):
    result_file = os.path.join(result_dir, os.path.basename(tf_bed))
    os.system("bedtools fisher -a %s -b %s -g %s > %s" % (mettl3_bed, tf_bed, chrom_size_file, result_file))


if __name__ == '__main__':
    pool = Pool()
    for tf_file in glob.glob("%s/*.bed" % tf_dir):
        pool.apply_async(fisher_test, (tf_file, ))
    pool.close()
    pool.join()
