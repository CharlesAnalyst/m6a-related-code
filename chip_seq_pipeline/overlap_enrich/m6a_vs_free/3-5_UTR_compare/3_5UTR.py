#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd
import subprocess
from multiprocessing import Pool

anno_type = "5UTR"    # 3UTR
result_dir = "/data5/galaxy/project/data/total_m6a_peak/annotate_peak/%s" % anno_type
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
species = "hg38"
peak_dir = "/data5/galaxy/project/data/total_m6a_peak"


def main():
    peak_list = glob.glob("%s/*.bed" % peak_dir)
    pool = Pool()
    for peak_file in peak_list:
        anno_result = os.path.join(result_dir, os.path.basename(peak_file).lower())
        pool.apply_async(get_annotated_region_peak, (peak_file, anno_result))
        # get_annotated_region_peak(peak_file, anno_result)
    pool.close()
    pool.join()


def get_annotated_region_peak(peak_bed, result_file):
    result_string = homer_annotate_peak(peak_bed)
    parse_annotated_region(result_string, result_file)


def homer_annotate_peak(peak_bed):
    command = "annotatePeaks.pl %s %s" % (peak_bed, species)
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    result_string = sub_p.communicate()[0]
    return result_string


def parse_annotated_region(result_string, result_file):
    line_list = str(result_string).strip().split("\n")
    col_names, content = line_list[0].split("\t"), line_list[1:]
    content_list = [line.split("\t") for line in content]
    df_content = pd.DataFrame(content_list, columns=col_names)
    # new start && annotation
    df_content["new_s"], df_content["new_a"] = (df_content["Start"].astype(int) - 1), (df_content["Annotation"].str.split("(").str[0].str.replace("' ", "")).str.strip()
    df_type = df_content.loc[df_content["new_a"] == anno_type]
    # print(df_type.head())
    df_bed = df_type[["Chr", "new_s", "End"]].sort_values(["Chr", "new_s"])
    df_bed.to_csv(result_file, sep="\t", header=False, index=False, quoting=False)


if __name__ == '__main__':
    main()