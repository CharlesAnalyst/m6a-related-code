#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd
from multiprocessing import Pool

############################################################
bed_dir = "/data5/galaxy/project/bam2bed_input/merged_data"
bed_files = glob.glob("%s/*.bed" % bed_dir)
peak_length = 374
result_dir = "/data5/galaxy/project/data/input_data/autosomal/input_bed/universal"
# complete_names = ["brain", "heart", "lung", "liver", "kidney", "placenta", "muscle", "stomach"]
############################################################


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def broke_peaks(i_bed):
    prefix = os.path.basename(i_bed).split(".bed")[0]
    # df = pd.read_table(i_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "a", "strand"])
    df = pd.read_table(i_bed, sep="\t", header=None, names=["chr", "start", "end"])
    df_peak = df[["chr", "start", "end"]]
    df_peak.loc[:, "length"] = df_peak["end"] - df_peak["start"]
    peak_list = []
    for name, values in df_peak.iterrows():
        group_list = list(chunks(range(values["start"], values["end"]), peak_length))
        for i_group_list in group_list:
            peak_list.append([values["chr"], i_group_list[0], (i_group_list[-1] + 1)])
    df_result = pd.DataFrame(peak_list)
    df_result.columns = ["chr", "start", "end"]
    result_bed = "%s/%s_discrete.bed" % (result_dir, prefix)
    df_result.to_csv(result_bed, sep="\t", header=None, index=False)


if __name__ == "__main__":
    broke_peaks("/data5/galaxy/project/data/input_data/autosomal/input_bed/universal/universal.bed")
    """
    pool = Pool(len(bed_files))
    for tissue_bed in bed_files:
        pool.apply_async(broke_peaks, (tissue_bed, ))
    pool.close()
    pool.join()
"""