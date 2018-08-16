#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import pandas as pd

#######################################################################
gtf_file = "/data4/database/GRCh38/GENCODE/gencode.v27.annotation.gff3"
result_dir = "/data5/galaxy/project/tf_analysis/TFBS_analysis/control"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = os.path.join(result_dir, "back_ground.bed")
###############################################################################


df = pd.read_table(gtf_file, comment="#", header=None)
df_gene = df[df.iloc[:, 2] == "gene"]
result_list = []
# consider positive and negative strand
for name, values in df_gene.iterrows():
    if values[6] == "+":
        if values[3] >= 2000:
            result_line = [values[0], (values[3] - 2000), values[4]]
        else:
            result_line = [values[0], values[3], values[4]]
    else:
        result_line = [values[0], values[3], (values[4] + 2000)]
    result_list.append(result_line)


df_result = pd.DataFrame(result_list)
df_result.to_csv(result_file, sep="\t", header=None, index=None)
