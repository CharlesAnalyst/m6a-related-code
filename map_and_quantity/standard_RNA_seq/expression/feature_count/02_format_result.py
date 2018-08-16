#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import pandas as pd


file = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/expression/featureCount/gene_counts.txt"
result_count_file = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/expression/featureCount/format_counts.txt"
result_length_file = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/expression/featureCount/length.txt"
result_col_file = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/expression/featureCount/colData.txt"

# only keep count and length
df = pd.read_table(file, sep="\t", comment="#", index_col=0)
old_cols = [x for x in df.columns if "bam" in x]
df_count = df.copy()[old_cols]
col_names = [x.split("/")[-1].split("_")[0].replace("-", "_") for x in df_count.columns]
df_count.columns = col_names
# df_count["Length"] = df["Length"]
df_count.to_csv(result_count_file, sep="\t", header=True, index=True)
#
df_len = pd.DataFrame(df["Length"])
df_len.to_csv(result_length_file, sep="\t", header=True, index=True)
#
condition_list = [x.split("_")[-1].lower() for x in col_names]
type_list = ["paired_end" for i in range(len(condition_list))]
df_col = pd.DataFrame({"condition": condition_list, "type": type_list})
df_col.index = col_names
df_col.to_csv(result_col_file, sep="\t", header=True, index=True)