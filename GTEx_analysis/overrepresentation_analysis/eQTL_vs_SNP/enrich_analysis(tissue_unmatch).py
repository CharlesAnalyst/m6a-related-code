#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import numpy as np
import glob
from scipy import stats
from multiprocessing import Pool

###############################################################################################################
eQTL_dir = "/data5/galaxy/project/snp_analysis/snp_db/GTEx_Analysis_v7_eQTL/format_bed/eQTL_snp/autosomal_snp/"
# eQTL_bed_1 = "%s/Heart_Left_Ventricle.bed" % eQTL_dir
# eQTL_bed_2 = "%s/Heart_Atrial_Appendage.bed" % eQTL_dir
# eQTL_list = [eQTL_bed_1, eQTL_bed_2]
eQTL_list = glob.glob("%s/*.bed" % eQTL_dir)
snp_bed = "/data/database/snp/autosomal_GRCh38_snp.bed"
#########################
m6a_dir = "/data5/galaxy/project/data/total_m6a_peak"
control_dir = "/data5/galaxy/project/data/input_data/input_cycle"
#########################
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_vs_snp/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = os.path.join(result_dir, "eQTL_vs_snp_no-match_result.txt")
#################################################################################################################


def get_control_dict():
    control_dict = {}
    control_list = glob.glob("%s/control-*_*.bed" % control_dir)
    for control in control_list:
        tissue = os.path.basename(control).split("_")[0].split("-")[1].lower()
        control_dict[tissue] = control_dict.get(tissue, []) + [control]
    return control_dict


def fisher_test(a, b, eQTL_bed, control_bed):
    c = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (eQTL_bed, control_bed)).read()
    d = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (snp_bed, control_bed)).read()
    oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
    return [oddsratio, pvalue]


def intersect_bed(eQTL_bed, m6a_bed, control_list):
    #          eQTL snp  HapMap
    # m6a          a       b
    # control      c       d
    a = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (eQTL_bed, m6a_bed)).read()
    b = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (snp_bed, m6a_bed)).read()
    # input bed
    pool, results = Pool(), []
    for control_bed in control_list:
        # [oddsratio, pvalue] = fisher_test(a, b, eQTL_bed, control_bed)
        result = pool.apply_async(fisher_test, (a, b, eQTL_bed, control_bed))
        results.append(result)
    result_list = [x.get() for x in results]
    return result_list


def main_method():
    control_dict, m6a_list = get_control_dict(), glob.glob("%s/*.bed" % m6a_dir)
    with open(result_file, 'w') as fw:
        fw.write("sub_tissue\tm6a_tissue\tfinal_or\n")
        for eQTL_bed in eQTL_list:
            sub_tissue = os.path.basename(eQTL_bed).split(".bed")[0]
            tissue = sub_tissue.split("_")[0].lower()
            print(tissue)
            for m6a_bed in m6a_list:
                m6a_tissue = os.path.basename(m6a_bed).split(".bed")[0].lower()
                control_list = control_dict[m6a_tissue]
                result_list = intersect_bed(eQTL_bed, m6a_bed, control_list)
                # filter
                oddsratio_list = [0.0 if float(line[1]) > 0.05 else line[0] for line in result_list]
                final_or = np.average(oddsratio_list)
                fw.write("%s\t%s\t%f\n" % (sub_tissue, m6a_tissue, final_or))


if __name__ == '__main__':
    main_method()
