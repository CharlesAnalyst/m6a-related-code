#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import pandas as pd
#

result_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/gene_bed"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
##########
genome_region_bed = "/data/database/GRCh38/GENCODE/Genes_ensembl.bed"
df_bed = pd.read_table(genome_region_bed, sep="\t", header=None, names=["chr", "a", "b", "gene", "c", "d"])
#####################
expre_file = "/data3/xs/tissue_m6a/2018.1/fig2/fig2_1_21/Total-nofilter.txt"
df = pd.read_table(expre_file, sep="\t", index_col=0)
df_exp = df[[col for col in df.columns if "exp" in col]]
del df_exp["expression_cv"]
for col in df_exp.columns:
    tissue = col.split("_exp")[0].lower()
    print(tissue)
    result_file = os.path.join(result_dir, "%s.bed" % tissue)
    df_col = df_exp[col]
    df = pd.DataFrame(df_col[df_col > 0]).reset_index()
    df_overlap = pd.merge(df, df_bed, how="left")
    df_bed = df_overlap[["chr", "a", "b", "gene", "c", "d"]].dropna(how="any").drop_duplicates()
    df_bed["a"], df_bed["b"] = df_bed["a"].astype(int), df_bed["b"].astype(int)
    df_sort = df_bed.sort_values(["chr", "a"])
    df_sort.to_csv(result_file, sep="\t", header=None, index=False)


