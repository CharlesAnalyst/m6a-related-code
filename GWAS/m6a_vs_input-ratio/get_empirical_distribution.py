#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
# function: generate empirical null distribution; re-generate random 100 control region;
# generate 100 randomly drawn size matched sets of control regions and determined the number of
# overlapping SNPs and determined the empirical distribution function of the background SNP overlap.


import os
import glob
import pandas as pd


#######################################################################################################################
intersect_dir = "/data5/galaxy/project/GWAS_analysis/cluster_by_map/ratio_version/intersect_snp_class/control"
result_dir = "/data5/galaxy/project/GWAS_analysis/cluster_by_map/ratio_version/empirical_distribution"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
count_result_file = os.path.join(result_dir, "empirical_distribution.txt")
#######################################################################################


def main():
    df = calculate_background_distribution()
    return df


def statistic_overlap_num(intersect_file):
    with open(intersect_file, 'r') as f:
        overlap_num = len(f.readlines())
    return overlap_num


def calculate_background_distribution():
    df, tissue_dict = pd.DataFrame(), {}
    bed_list = glob.glob("%s/*.bed" % intersect_dir)
    for i_bed in bed_list:
        snp_trait = "_".join(os.path.basename(i_bed).split("_")[:-1])
        tissue_dict[snp_trait] = tissue_dict.get(snp_trait, []) + [i_bed]
    for tissue, bed_list in tissue_dict.items():
        # print(tissue)
        num_list = []
        for i_bed in bed_list:
            i_num = statistic_overlap_num(i_bed)
            num_list.append(i_num)
        df[tissue] = num_list
    # df.to_csv(count_result_file, sep="\t", index=False)
    return df


# if __name__ == "__main__":
#     calculate_background_distribution()
