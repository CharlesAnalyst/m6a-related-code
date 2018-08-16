#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
import permutation_test as pt

# compare each gene by each tissue
# no zscore; just filter by cutoff and RPKM
####################################################################################################################################
fimo_table = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/methylation_complex_vs_others/filter_by_cutoff-RPKM.txt"
id_transform_file = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/04_intersect/id_trans.txt"
complex_dir = "/data5/galaxy/project/what/TF_narrowPeak/clean_data"
complex_gene_list = [os.path.basename(x).split(".bed")[0].upper() for x in glob.glob("%s/*.bed" % complex_dir)]
result_file = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/methylation_complex_vs_others/u-test_result.txt"
#############################################################################################################

#
trans_id = pd.read_table(id_transform_file, sep="\t")
df_motif = pd.read_table(fimo_table, sep="\t")
df_merge = pd.merge(df_motif, trans_id, left_on="ensembl_id", right_on="ENSEMBL", how="left")
# df_merge = df_merge.dropna(how="any")
del df_merge["ENSEMBL"]
del df_merge["ensembl_id"]
df = df_merge.set_index(["SYMBOL"])
print(len(df))

with open(result_file, 'w') as fw:
    fw.write("tissue\tstatistic\tpvalue\ttarget_mean\tcontrol_mean\n")
    for tissue in df.columns:
        print(tissue)
        df_sub = pd.DataFrame(df[tissue])
        target_gene_list, control_gene_list = [], []
        for name, value in df_sub.iterrows():
            if float(value) > 0:   # remove zero
                if name in complex_gene_list:
                    target_gene_list.append(float(value))
                else:
                    control_gene_list.append(float(value))
        print(target_gene_list)
        statistic, pvalue = mannwhitneyu(target_gene_list, control_gene_list)
        # pvalue = pt.permutationtest(target_gene_list, control_gene_list)
        # statistic, pvalue = ttest_ind(target_gene_list, control_gene_list)
        target_aver, control_aver = np.mean(target_gene_list), np.mean(control_gene_list)
        fw.write("%s\t%f\t%f\t%s\t%s\n" % (tissue, statistic, pvalue, target_aver, control_aver))
        if tissue == "placenta":
            print(target_gene_list)
            print("!!!!!!!!!!")
            print(control_gene_list)
        # fw.write("%s\t%f\n" % (tissue, pvalue))
