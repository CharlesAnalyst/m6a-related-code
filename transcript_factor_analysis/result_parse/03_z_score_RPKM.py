#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import math
import numpy
import pandas as pd
from pandas import DataFrame


def main_method():
    work_dir = "/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result"
    table = "%s/bigTable.txt" % work_dir
    z_score_table = "%s/bigTable_zscore.txt" % work_dir

    with open(z_score_table, 'w') as fw:
        with open(table, 'r') as f:
            title = f.readline()
            fw.write(title)
            for line in f.readlines():
                motif_id, values = line.split("\t")[0], [float(value) for value in line.split("\t")[1:]]
                df = pd.DataFrame(values, columns=["ratio"])
                new_list = list((df["ratio"] - df["ratio"].mean())/df["ratio"].std())
                score_list = [str(i) for i in new_list]
                # To correct for samples not enriched over background,
                # we set the z-scores for motifs with an enrichment ratio below 1.2 to 0 in that particular sample.
                # ratio_less_cutoff_index = [values.index(ratio) for ratio in values if ratio < 1.2]
                # for i in ratio_less_cutoff_index:
                #    score_list[i] = "0.0"
                new_line = "%s\t%s" % (motif_id, "\t".join(score_list))
                fw.write(new_line.strip() + "\n")


main_method()
