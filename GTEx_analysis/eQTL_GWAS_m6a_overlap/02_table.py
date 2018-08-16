#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import random
import subprocess
import numpy as np
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


#######
epoch = 1000
######
total_eQTL = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/format_bed/total_eQTL_snp/total_eQTL_snp.bed"
df = pd.read_table(total_eQTL, sep="\t", header=None, names=["a", "b", "c", "rs_id"])
total_eQTL_list = list(set(df["rs_id"].tolist()))
######
gwas_snp = "/data5/galaxy/project/GWAS_db/03_replicated_GRCh38/all_replicated_gwas_snp.bed"
df = pd.read_table(gwas_snp, sep="\t", header=None, names=["a", "b", "c", "rs_id"])
gwas_string = str(df.to_string(header=False, index=False, col_space=4))
gwas_string = "\n".join(["\t".join(x.split()) for x in gwas_string.split("\n")])
######
input_bed = "/data5/galaxy/project/data/input_data/autosomal/universal/input_bed/universal.bed"
command = "bedtools intersect -a stdin -b %s -wa" % input_bed
sub_pro = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
gwas_line_list = sub_pro.communicate(gwas_string)[0].strip().split("\n")
# total_gwas_list = list(set(df["rs_id"].tolist()))
total_gwas_list = [x.split("\t")[3] for x in gwas_line_list]
###############
overlap_set = set(total_eQTL_list).intersection(set(total_gwas_list))
print(len(overlap_set))
#


def plot_histogram(data_list):
    num_bins, mean, std = 20, np.mean(data_list), np.std(data_list)
    n, bins, patches = plt.hist(data_list, num_bins, normed=1, facecolor="blue", alpha=0.5)
    y = mlab.normpdf(bins, mean, std)
    plt.plot(bins, y, "r--")
    plt.xlabel("Number of overlapped snp")
    plt.ylabel("Probability")
    plt.title(r"Histogram of IQ: $mean=%s$, $std=%s$" % (str(mean), str(std)))
    plt.subplots_adjust(left=0.15)
    plt.show()


len_list = []
for i in range(epoch):
    selected_eQTL_list = random.sample(total_eQTL_list, 24684)
    # print(len(selected_eQTL_list))
    selected_gwas_list = random.sample(total_gwas_list, 379)
    # print(len(selected_gwas_list))
    overlap_set = set(selected_eQTL_list).intersection(set(selected_gwas_list))
    # print(overlap_set)
    len_list.append(len(overlap_set))
plot_histogram(len_list)