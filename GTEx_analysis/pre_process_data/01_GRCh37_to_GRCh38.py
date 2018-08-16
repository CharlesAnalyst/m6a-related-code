#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob

#####################################################################################################################
in_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/01_2_format_bed_GRCh37"
postfix = ".bed"
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/02_bed_GRCh38"
chain_file = "/data/database/liftover_chain/hg19ToHg38.over.chain"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#######################################################################################################################


in_list = glob.glob("%s/*%s" % (in_dir, postfix))
for in_bed in in_list:
    out_bed = os.path.join(result_dir, os.path.basename(in_bed))
    un_lift = os.path.join(result_dir, "un_lift.bed")
    os.system("liftOver %s %s %s %s" % (in_bed, chain_file, out_bed, un_lift))
    os.system("rm %s" % un_lift)
