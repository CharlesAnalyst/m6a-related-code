#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd
import subprocess
from multiprocessing import Pool

species = "hg38"
peak_dir = "/data5/galaxy/project/methyl_m6a/data/diff_m6a_peak/macs2_bdgdiff/permutation/"
# result_dir = "/data5/galaxy/project/methyl_m6a/data/diff_m6a_peak/macs2_bdgdiff/permutation/"


def main():
    peak_list = glob.glob("%s/*/*.bed" % peak_dir)
    # peak_file = "%s/DTK_vs_CGR8_c3.0_cond1.bed" % peak_dir
    pool = Pool(processes=20)
    for peak_file in peak_list:
        pool.apply_async(single_run, (peak_file, ))
    pool.close()
    pool.join()
    # for peak_file in peak_list:
    #     single_run(peak_file)


def single_run(peak_file):
    print(os.path.basename(peak_file))
    result_string = homer_annotate_peak(peak_file)
    # anno_result = os.path.join(result_dir, os.path.basename(peak_file).replace(".bed", "_anno.txt"))
    anno_result = peak_file.replace(".bed", "_anno.txt")
    # genes_ensemble = get_annotated_genes(result_string)
    df = get_annotated_genes(result_string)
    write_to_file(df, anno_result)


def homer_annotate_peak(peak_bed):
    command = "annotatePeaks.pl %s %s" % (peak_bed, species)
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    result_string = sub_p.communicate()[0]
    return result_string


def get_annotated_genes(result_string):
    line_list = str(result_string).strip().split("\\n")
    # print(len(line_list))
    col_names, content = line_list[0].split("\\t"), line_list[1:]
    # print(col_names)
    content_list = [line.split("\\t") for line in content]
    df_content = pd.DataFrame(content_list, columns=col_names)
    # print(df_content.head())
    # genes_ensemble = df_content["Nearest Ensembl"]
    df = df_content[["Nearest Ensembl", "Peak Score", "Annotation"]]  # .replace(r"\s+", np.nan, regex=True).dropna()
    df = df[(df["Annotation"] != "Intergenic") & (df["Nearest Ensembl"] != "")]
    df["Peak Score"] = df["Peak Score"].astype(float)
    df = df[["Nearest Ensembl", "Peak Score"]].sort_values(["Nearest Ensembl"])
    # print(df.head(n=30))
    return df


def write_to_file(df, result_file):
    # with open(result_file, 'w') as fw:
    #     fw.writelines(["%s\n" % x.strip() for x in result_list if len(x) > 0])
    df.to_csv(result_file, sep="\t", header=False, index=False)


if __name__ == '__main__':
    main()