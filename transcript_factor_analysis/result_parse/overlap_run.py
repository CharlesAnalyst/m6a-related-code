#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool


####################################################################################################
tfbs_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_hg19_to_GRCh38"
cycle_dir = "/data5/galaxy/project/data/input_data/input_cycle"
m6a_dir = "/data5/galaxy/project/data/total_m6a_peak"
intersect_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/04_intersect"
if not os.path.exists(intersect_dir):
    os.makedirs(intersect_dir)
####################################################################################################


def overlap_input(input_bed, tfbs_bed, tissue_name):
    result_dir = "%s/%s/control" % (intersect_dir, tissue_name)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    cyc_num = os.path.basename(input_bed).split("_")[1].split(".bed")[0]
    result_file = os.path.join(result_dir, "%s-%s" % (cyc_num, os.path.basename(tfbs_bed)))
    os.system("bedtools intersect -a %s -b %s > %s" % (input_bed, tfbs_bed, result_file))


def overlap_ip(m6a_bed, tfbs_bed, tissue_name):
    result_dir = "%s/%s/ip" % (intersect_dir, tissue_name)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    result_file = os.path.join(result_dir, os.path.basename(tfbs_bed))
    os.system("bedtools intersect -a %s -b %s > %s" % (m6a_bed, tfbs_bed, result_file))


def get_bed_according_tfbs(name, in_list):
    out_file = ""
    for x in in_list:
        if name in x:
            out_file = x
            break
    return out_file


def main_method():
    tfbs_list = glob.glob("%s/*.bed" % tfbs_dir)
    m6a_list = glob.glob("%s/*.bed" % m6a_dir)
    for tfbs in tfbs_list:
        print(tfbs)
        for m6a_bed in m6a_list:
            print(m6a_bed)
            tissue_name = os.path.basename(m6a_bed).split(".bed")[0].lower()
            input_list = glob.glob("%s/control-%s_*.bed" % (cycle_dir, tissue_name))
            print(os.path.basename(tfbs), os.path.basename(m6a_bed), input_list[0])
            overlap_ip(m6a_bed, tfbs, tissue_name)
            pool = Pool()
            for input_bed in input_list:
                # overlap_input(input_bed, m6a_bed)
                pool.apply_async(overlap_input, (input_bed, tfbs, tissue_name))
            pool.close()
            pool.join()


if __name__ == "__main__":
    main_method()
