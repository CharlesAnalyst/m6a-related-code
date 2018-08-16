#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import pandas as pd

####################################################################################
bed_file = "/data5/galaxy/project/tf_analysis/TFBS_analysis/control/back_ground.bed"
peak_length = 504  # [381.31795456538902, 429.01515375517448, 586.99684030558137, 507.60837378640775, 618.1096905947519]
result_file = "/data5/galaxy/project/tf_analysis/TFBS_analysis/control/bg_discrete.bed"
########################################################################################


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def broke_peaks(i_bed):
    df = pd.read_table(i_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "a", "strand"])
    df_peak = df[["chr", "start", "end"]]
    df_peak.loc[:, "length"] = df_peak["end"] - df_peak["start"]
    peak_list = []
    for name, values in df_peak.iterrows():
        group_list = list(chunks(range(values["start"], values["end"]), peak_length))
        for i_group_list in group_list:
            peak_list.append([values["chr"], i_group_list[0], (i_group_list[-1] + 1)])
    df_result = pd.DataFrame(peak_list)
    df_result.columns = ["chr", "start", "end"]
    df_result.to_csv(result_file, sep="\t", header=None, index=False)


if __name__ == "__main__":
    broke_peaks(bed_file)
