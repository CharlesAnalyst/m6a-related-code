#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool


################
base_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/macs2_peak"
result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/diff_peak"
###
compare_list = ["DTK:CGR8"]
tissue_prefix_list = [os.path.basename(x).split("_peaks.xls")[0] for x in glob.glob("%s/*_peaks.xls" % base_dir)]
################


def process_each_pair(pair_list):
    tissue_1, tissue_2 = pair_list[0], pair_list[1]
    print(tissue_1, tissue_2)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    treat_bdg_1, control_bdg_1, final_num_1 = get_effective_seq_depth(tissue_1)
    treat_bdg_2, control_bdg_2, final_num_2 = get_effective_seq_depth(tissue_2)
    os.system("macs2 bdgdiff --t1 %s --c1 %s --t2 %s --c2 %s --d1 %f --d2 %f --outdir %s --o-prefix '%s_vs_%s'" %
              (treat_bdg_1, control_bdg_1, treat_bdg_2, control_bdg_2, final_num_1, final_num_2, result_dir, tissue_1, tissue_2))


def get_effective_seq_depth(tissue_prefix):
    treat_bdg = "%s/%s_treat_pileup.bdg" % (base_dir, tissue_prefix)
    control_bdg = "%s/%s_control_lambda.bdg" % (base_dir, tissue_prefix)
    xls = "%s/%s_peaks.xls" % (base_dir, tissue_prefix)
    print(xls)
    return_str = os.popen('egrep "tags after filtering in treatment|tags after filtering in control" %s' % xls)
    line_list, number_list = return_str.read().strip().split("\n"), []
    for line in line_list:
        number_list.append(float(line.split(":")[1].strip()))
    final_num = min(number_list)
    return treat_bdg, control_bdg, final_num


if __name__ == "__main__":
    pool = Pool(processes=5)
    for x in compare_list:
        paired_tissue = x.split(":")
        # process_each_pair(paired_tissue)
        pool.apply_async(process_each_pair, (paired_tissue,))
    pool.close()
    pool.join()
