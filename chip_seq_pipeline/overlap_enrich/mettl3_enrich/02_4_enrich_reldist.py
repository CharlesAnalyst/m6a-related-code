#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
from multiprocessing import Pool


tf_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean"
mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
result_dir = "/data5/galaxy/project/mettl3_enrich/result_parse/bedtools_reldist"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def reldist_test(tf_bed):
    result_file = os.path.join(result_dir, os.path.basename(tf_bed))
    os.system("bedtools reldist -a %s -b %s > %s" % (mettl3_bed, tf_bed, result_file))


if __name__ == '__main__':
    pool = Pool()
    for tf_file in glob.glob("%s/*.bed" % tf_dir):
        pool.apply_async(reldist_test, (tf_file, ))
    pool.close()
    pool.join()
