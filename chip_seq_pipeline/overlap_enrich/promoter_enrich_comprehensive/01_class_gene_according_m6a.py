#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import subprocess
import pandas as pd
from multiprocessing import Pool


###########################################################################################
m6a_bed_dir = "/data5/galaxy/project/data/total_m6a_peak"
m6a_list = glob.glob("%s/*.bed" % m6a_bed_dir)
#
original_gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/gene_bed"
gene_list = glob.glob("%s/*.bed" % original_gene_dir)
original_promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/promoter_bed"
promoter_list = glob.glob("%s/*.bed" % original_promoter_dir)
#
gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_gene/gene_bed"
if not os.path.exists(gene_dir):
    os.makedirs(gene_dir)
###########################################################################################


def process_by_each_tissue(gene_bed, gene_promoter_bed, m6a_bed):
    df_m6a_genes = get_gene_according_m6a(gene_bed, gene_promoter_bed, m6a_bed)


def get_gene_according_m6a(gene_bed, gene_promoter_bed, m6a_bed):
    gene_bed, gene_promoter_bed = get_gene_promoter(gene_bed, gene_promoter_bed)
    df_m6a_genes = get_m6a_gene(gene_bed, gene_promoter_bed, m6a_bed)
    return df_m6a_genes


def get_gene_promoter(gene_bed, gene_promoter_bed):
    df = pd.read_table(gene_bed, sep="\t", header=None, names=["a", "b", "c", "name", "d", "e"])
    df_pro = pd.read_table(gene_promoter_bed, sep="\t", header=None, names=["f", "g", "h", "name", "i", "j"])
    df_merge = pd.merge(df, df_pro, on="name", how="right").sort_values(["a", "b"])
    df_merge["start"], df_merge["end"] = df_merge.min(axis=1, numeric_only=True), df_merge.max(axis=1, numeric_only=True)
    df_gene_pro = df_merge[["a", "start", "end", "name", "d", "e"]]
    print(df_gene_pro.head())
    return df, df_gene_pro


def get_m6a_gene(df_gene, df_gene_pro, m6a_bed):
    gene_pro_string = str(df_gene_pro.to_string(header=False, index=False, col_space=4))
    gene_pro_string = "\n".join(["\t".join(x.split()) for x in gene_pro_string.split("\n")])
    m6a_gene_bed = os.path.join(gene_dir, os.path.basename(m6a_bed))
    command = "bedtools intersect -a stdin -b %s -wa | sort -k1,1 -k2,2n | uniq" % m6a_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    m6a_gene_promoter = str(sub_p.communicate(gene_pro_string)[0]).strip().split("\n")
    # exclude promoter
    gene_names = []
    for line in m6a_gene_promoter:
        gene_names.append(line.split("\t")[3])
    df_m6a_genes = df_gene[df_gene["name"].isin(gene_names)]
    df_m6a_genes.to_csv(m6a_gene_bed, sep="\t", header=None, index=False)
    return df_m6a_genes


if __name__ == '__main__':
    gene_bed_dict, gene_promoter_bed_dict = {}, {}
    for g_bed in gene_list:
        gene_bed_dict[os.path.basename(g_bed).split(".bed")[0].lower()] = g_bed
    for g_p_bed in promoter_list:
        gene_promoter_bed_dict[os.path.basename(g_p_bed).split(".bed")[0].lower()] = g_p_bed
    pool = Pool()
    for i_bed in m6a_list:
        tissue = os.path.basename(i_bed).split(".bed")[0].lower()
        g_bed = gene_bed_dict[tissue]
        g_p_bed = gene_bed_dict[tissue]
        pool.apply_async(process_by_each_tissue, (g_bed, g_p_bed, i_bed))
        # process_by_each_tissue(g_bed, i_bed)
    pool.close()
    pool.join()
    print("all done！")
