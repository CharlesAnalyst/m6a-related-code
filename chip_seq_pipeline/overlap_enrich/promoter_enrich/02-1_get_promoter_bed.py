#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd
from multiprocessing import Pool


########
gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_gene/gene_bed"
promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/m6a_gene/promoter_bed"
promoter_bed = "/data5/galaxy/project/data/promoter/2k_100/genes_promoters.bed"
###################


gene_list = glob.glob("%s/*.bed" % gene_dir)
gene_bed_dict = {}
for gene_bed in gene_list:
    tissue_name = os.path.basename(gene_bed).split(".bed")[0].lower()
    gene_bed_dict[tissue_name] = gene_bed

df_promoter = pd.read_table(promoter_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand"])
for tissue_name, gene_bed in gene_bed_dict.items():
    print(tissue_name)
    result_bed = os.path.join(promoter_dir, os.path.basename(gene_bed))
    df_gene = pd.read_table(gene_bed, sep="\t", header=None, names=["a", "b", "c", "name", "e", "f"])
    df_merge = pd.merge(df_gene, df_promoter, on="name", how="left")
    for col in ["a", "b", "c", "e", "f"]:
        del df_merge[col]
    df_bed = df_merge.dropna(how="any").drop_duplicates()
    df_bed["start"] = df_bed["start"].astype(int)
    df_bed["end"] = df_bed["end"].astype(int)
    df_bed["score"] = df_bed["score"].astype(int)
    df_result = df_bed[["chr", "start", "end", "name", "score", "strand"]]
    df_result.to_csv(result_bed, sep="\t", header=False, index=False, quoting=False)