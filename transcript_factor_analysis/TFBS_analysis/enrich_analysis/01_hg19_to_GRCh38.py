#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool

bed_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/02_combined_tf"
result_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_hg19_to_GRCh38"


def liftover_version_hg19(in_file):
    prefix = os.path.basename(in_file).split(".bed")[0].upper()
    hg38_bed, unmap_bed = os.path.join(result_dir, "hg38-%s.bed" % prefix), os.path.join(result_dir,
                                                                                         "unmap-%s.bed" % prefix)
    os.system("liftOver %s /home/galaxy/software/liftover/hg19ToHg38.over.chain %s %s"
              % (in_file, hg38_bed, unmap_bed))
    os.remove(unmap_bed)
    return hg38_bed


def pick_chr1_22_XY(hg38_bed):
    result_file = os.path.join(result_dir, os.path.basename(hg38_bed).split("hg38-")[1])
    os.system("less %s | awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y)$/' > %s" % (hg38_bed, result_file))
    os.remove(hg38_bed)


def transform_data(bed):
    lift_bed = liftover_version_hg19(bed)
    pick_chr1_22_XY(lift_bed)


if __name__ == '__main__':
    pool = Pool()
    pool.map_async(transform_data, glob.glob("%s/*.bed" % bed_dir))
    pool.close()
    pool.join()







