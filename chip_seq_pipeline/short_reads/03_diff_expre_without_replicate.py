#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import pandas as pd

########################################################################################################################
count_file = "/data5/galaxy/project/DNMT1_KO/featureCounts/rpkm_input.csv"
result_dir = "/data5/galaxy/project/DNMT1_KO/diff_expre"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
PROPORTATION = 0.1
result_up_gene, result_down_gene = os.path.join(result_dir, "up_genes.txt"), os.path.join(result_dir, "down_genes.txt")
#################################################################


df = pd.read_csv(count_file, index_col=0)
df["diff_rpkm"] = df.loc[:, "DTK_In"] - df.loc[:, "CGR8_In"]
df_up_raw, df_down_raw = df[df["diff_rpkm"] > 0], df[df["diff_rpkm"] <= 0]  # 10583 # 41967
print(df_up_raw.head())
up_len, down_len = len(df_up_raw), len(df_down_raw)
print(up_len, down_len)
query_up_len, query_down_len = int(up_len * PROPORTATION), int(down_len * PROPORTATION)
print(query_up_len, query_down_len)
df_up = df_up_raw.sort_values(["diff_rpkm"], ascending=False).iloc[:query_up_len, :]
print(df_up.head())
df_down = df_down_raw.sort_values(["diff_rpkm"], ascending=True).iloc[:query_down_len, :]
up_genes, down_gens = df_up.index, df_down.index
with open(result_up_gene, 'w') as fw:
    fw.writelines(["%s\n" % gene for gene in up_genes])
with open(result_down_gene, 'w') as fw:
    fw.writelines(["%s\n" % gene for gene in down_gens])