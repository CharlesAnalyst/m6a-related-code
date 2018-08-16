#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
from collections import Counter


#############################################################################
fimo_dir = "/data5/galaxy/project/tf_pathway/fimo_result"
intersect_snp_dir = "/data5/galaxy/project/tf_pathway/m6a_snp_intersect/snp/"
snp_list = glob.glob("%s/*.bed" % intersect_snp_dir)
result_dir = "/data5/galaxy/project/tf_pathway/m6aSNP_TF_intersect"
result_tf_dir, result_snp_dir, tf_snp_dir = "%s/TF" % result_dir, "%s/snp" % result_dir, "%s/tf_snp" % result_dir
for i_dir in [result_dir, result_tf_dir, result_snp_dir, tf_snp_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)
###################################################


def process_single_tissue(tissue):
    # os.chdir(os.path.join(fimo_dir, tissue.capitalize()))
    df = pd.read_table("/data5/galaxy/project/tf_pathway/fimo_result/all_tissues/fimo.txt", sep="\t")
    df["seq_start"] = df["sequence_name"].str.split(":").str[1].str.split("-").str[0]
    df["motif_start"] = df["seq_start"].astype(int) + df["start"].astype(int)
    df["motif_stop"] = df["seq_start"].astype(int) + df["stop"].astype(int) + 1
    df["chr"] = df["sequence_name"].str.split(":").str[0]
    df["motif_name"] = df["# motif_id"].astype(str) + ":" + df["motif_start"].astype(str) + "-" + df["motif_stop"].astype(str)
    df_bed = df[["chr", "motif_start", "motif_stop", "motif_name"]]
    # df_negative = df[df.strand == "-"]
    df_uniq = df_bed.drop_duplicates()
    df_uniq.to_csv("motifs.bed", sep="\t", header=None, index=False)
    intersect_result = "%s/raw_%s.txt" % (result_dir, tissue)
    snp_bed = get_intersect_snp_by_tissue(tissue)
    # snp_bed = ""
    print(snp_bed)
    os.system("bedtools intersect -a %s -b motifs.bed -wa -wb > %s" % (snp_bed, intersect_result))
    get_intersect_tf_list(tissue, intersect_result)
    get_intersect_snp_list(tissue, intersect_result)
    # return intersect_result


def get_intersect_tf_list(tissue, raw_file):
    df = pd.read_table(raw_file, sep="\t", header=None)
    df["tf"] = df.iloc[:, 7].str.split(":").str[0]
    ########
    tf_snp_file = "%s/%s.txt" % (tf_snp_dir, tissue)
    with open(tf_snp_file, 'w') as fw:
        for name, values in df.iterrows():
            fw.write("%s\t%s\n" % (values[7].split(":")[0], values[3]))
    ########
    tf_list = df["tf"].tolist()
    # count snp number overlapped with each transcript factor
    count = Counter(tf_list)
    sorted_count = sorted(count.items(), key=lambda a: a[1], reverse=True)
    result_file = "%s/%s.txt" % (result_tf_dir, tissue)
    with open(result_file, 'w') as fw:
        fw.write("transcript factor\tsnp number\n")
        for x, c in sorted_count:
            fw.write(x + "\t" + str(c) + "\n")


def get_intersect_snp_list(tissue, raw_file):
    df = pd.read_table(raw_file, sep="\t", header=None)
    df_snp_bed = df.iloc[:, 0:4].drop_duplicates()
    result_file = "%s/%s.bed" % (result_snp_dir, tissue)
    df_snp_bed.to_csv(result_file, sep="\t", header=None, index=False)


def get_intersect_snp_by_tissue(tissue):
    result_snp = ""
    for snp_bed in snp_list:
        if tissue in snp_bed.lower():
            result_snp = snp_bed
            break
    return result_snp


def process_all():
    # tissue_list = [x.split("/")[-2].lower() for x in glob.glob("%s/*/fimo.txt" % fimo_dir)]
    # for t in tissue_list:
    #    print(t)
    process_single_tissue("all_tissues")


if __name__ == '__main__':
    process_all()
