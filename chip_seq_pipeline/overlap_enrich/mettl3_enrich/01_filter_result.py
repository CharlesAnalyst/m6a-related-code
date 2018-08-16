#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import pandas as pd

in_file = "/data5/galaxy/project/mettl3_enrich/result_parse/overlap-statistical-test_result.txt"
result_file = "/data5/galaxy/project/mettl3_enrich/result_parse/sig_result.txt"
df = pd.read_table(in_file, sep="\t")
df_sig = df[df["p.adjust"] <= 0.05]
df_sig.to_csv(result_file, sep="\t", index=False, quoting=False)

