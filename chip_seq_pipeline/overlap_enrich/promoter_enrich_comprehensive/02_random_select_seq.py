#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import random
import pandas as pd

######################
subsample_size = 100
epoch = 10
total_promoter_bed = "/data5/galaxy/project/data/promoter/human/2k_100/genes_promoters.bed"
##############################################################################################################
base_dir = "/data5/galaxy/project/promoter_TF_enrich/data"
sub_dir = "m6a_gene"  # "total_gene"
gene_dir = "%s/%s/gene_bed" % (base_dir, sub_dir)
promoter_dir = "%s/%s/promoter_bed" % (base_dir, sub_dir)
###
work_dir = "%s/%s" % (base_dir, "random_%s" % sub_dir.split("_gene")[0])
result_gene_dir, result_pro_dir = os.path.join(work_dir, "gene_bed"), os.path.join(work_dir, "promoter_bed")
for i_dir in [result_gene_dir, result_pro_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)
###############################################################


def main():
    gene_list = glob.glob("%s/*.bed" % gene_dir)
    for i_gene in gene_list:
        print(os.path.basename(i_gene))
        gene_list = filter_gene_by_promoter(i_gene)
        random_select_gene(i_gene, gene_list)


def filter_gene_by_promoter(gene_bed):
    promoter_bed = os.path.join(promoter_dir, os.path.basename(gene_bed))
    name_list = pd.read_table(promoter_bed, sep="\t", header=None).iloc[:, 3].tolist()
    with open(gene_bed, 'r') as f:
        gene_bed_list = f.readlines()
        gene_list = [line for line in gene_bed_list if line.split("\t")[3] in name_list]
    return gene_list


def random_select_gene(gene_bed, gene_list):
    for i in range(epoch):
        selected_gene_list = random.sample(gene_list, subsample_size)
        #
        gene_result = os.path.join(result_gene_dir, "%d-%s" % ((i + 1), os.path.basename(gene_bed)))
        get_gene(gene_result, selected_gene_list)
        #
        pro_result = os.path.join(result_pro_dir, "%d-%s" % ((i + 1), os.path.basename(gene_bed)))
        get_promoter(selected_gene_list, pro_result)


def get_gene(gene_result, selected_gene_list):
    with open(gene_result, 'w') as fw:
        fw.writelines(selected_gene_list)


def get_promoter(gene_list, promoter_result_bed):
    df_gene_bed = pd.DataFrame([line.strip().split() for line in gene_list])
    df_gene_bed.columns = ["a", "b", "c", "name", "e", "f"]
    df_promoter = pd.read_table(total_promoter_bed, sep="\t", header=None,
                                names=["chr", "start", "end", "name", "score", "strand"])
    df_merge = pd.merge(df_gene_bed, df_promoter, on="name", how="left")
    # print(df_merge.head())
    for col in ["a", "b", "c", "e", "f"]:
        del df_merge[col]
    df_bed = df_merge.dropna(how="any").drop_duplicates()
    df_bed["start"] = df_bed["start"].astype(int)
    df_bed["end"] = df_bed["end"].astype(int)
    df_bed["score"] = df_bed["score"].astype(int)
    df_promoter_bed = df_bed[["chr", "start", "end", "name", "score", "strand"]]
    df_promoter_bed.to_csv(promoter_result_bed, sep="\t", header=None, index=False)


if __name__ == '__main__':
    main()
