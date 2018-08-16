#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import pandas as pd


def main_method():
    work_dir = "/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result"
    raw_file = "%s/bigTable.txt" % work_dir
    current_file = "%s/bigTable_zscore.txt" % work_dir
    result_file = "%s/zscore_filter-by-cutoff.txt" % work_dir
    df_raw = pd.read_table(raw_file, sep="\t", index_col=0)
    df_cur = pd.read_table(current_file, sep="\t", index_col=0)

    with open(result_file, 'w') as fw:
        fw.write("motif_id\t%s\n" % "\t".join(df_raw.columns))
        for name, values in df_cur.iterrows():
            values_raw = df_raw.loc[name, :]
            line_list = []
            for i in range(len(df_raw.columns)):
                if values_raw[i] > 1.1065:
                    line_list.append(values[i])
                else:
                    line_list.append(0.0)
            fw.write("%s\t%s\n" % (name, "\t".join([str(x) for x in line_list])))


main_method()
