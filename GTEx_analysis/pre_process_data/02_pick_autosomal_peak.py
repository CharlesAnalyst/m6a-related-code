#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob

#################################################################################################################
eQTL_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/02_bed_GRCh38"
bed_list = glob.glob("%s/*.bed" % eQTL_dir)
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/02_auto_bed_GRCh38"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#################################################################################################################


def main():
    for x in bed_list:
        pick_chr1_22(x)


def pick_chr1_22(in_bed):
    result_file = os.path.join(result_dir, os.path.basename(in_bed))
    # os.system("less %s | awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y)$/' > %s" % (in_bed, result_file))
    os.system("less %s | awk '$1 ~ /^chr(1?[0-9]|2[0-2])$/' | sort -k1,1 -k2,2n | uniq > %s" % (in_bed, result_file))
    # os.remove(in_bed)


if __name__ == '__main__':
    main()