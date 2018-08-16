#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import numpy as np
import datetime
from scipy import stats
import get_empirical_distribution
import calculate_pvalue


###################################################################################################################
treat_dir = "/data5/galaxy/project/GWAS_analysis/cluster_by_map/ratio_version/intersect_snp_class/condition"
control_dir = "/data5/galaxy/project/GWAS_analysis/cluster_by_map/ratio_version/intersect_snp_class/control"
result_dir = "/data5/galaxy/project/GWAS_analysis/cluster_by_map/ratio_version/treat_vs_control_ratio"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = "%s%sratio_static.txt" % (result_dir, os.sep)
###################################################################################################################


def main():
    result_list, control_overlapped_num_list = [], []
    df = get_empirical_distribution.main()
    treat_list, control_list = glob.glob("%s/*.bed" % treat_dir), glob.glob("%s/*.bed" % control_dir)
    treat_dict, control_dict = generate_map_dict(treat_list, "treat"), generate_map_dict(control_list, "control")
    for snp_trait, bed_list in treat_dict.items():
        treat_overlapped_num = statistic_overlap_num(bed_list[0])
        try:
            control_overlapped_num_list = df[snp_trait].tolist()
        except KeyError:
            print("%s not in df(empirical distribution table)!" % snp_trait)
        medium_ratio = np.mean([((treat_overlapped_num+0.5) / (x+0.5)) for x in control_overlapped_num_list])
        pvalue = calculate_pvalue.calculate_pvalue(treat_overlapped_num, control_overlapped_num_list)
        print("%s\t%f\t%f" % (snp_trait, medium_ratio, pvalue))
        result_list.append("%s\t%f\t%f\n" % (snp_trait, medium_ratio, pvalue))
    write_to_file(result_list)


def generate_map_dict(bed_list, data_type):  # treat or control
    map_dict = {}
    for bed in bed_list:
        snp_trait = ""
        if data_type == "treat":    # inflammatory-measurement_testis.bed
            snp_trait = os.path.basename(bed).split(".bed")[0]
        elif data_type == "control":    # response-to-drug_stomach_97.bed
            snp_trait = "_".join(os.path.basename(bed).split("_")[:-1])
        map_dict[snp_trait] = map_dict.get(snp_trait, []) + [bed]
    return map_dict


def statistic_overlap_num(intersect_file):
    with open(intersect_file, 'r') as f:
        overlap_num = len(f.readlines())
    return overlap_num


def write_to_file(result_list):
    with open(result_file, 'w') as fw:
        fw.write("snp_trait\tmedium_ratio\tpvalue\n")
        fw.writelines(result_list)


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(start_time)
    #
    main()
    #
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
