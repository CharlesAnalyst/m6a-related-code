#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import re
import math
import subprocess
import pandas as pd
import numpy as np
from scipy import stats


#############################################################################################
expression_file = "/data5/galaxy/project/DNMT1_KO/featureCounts/rpkm_input.csv"
up_peak_gene = "/data5/galaxy/project/DNMT1_KO/diff_peak/DTK_vs_CGR8_c3.0_cond1_anno.txt"
reference_genome = "/data/database/GRCm38/GENCODE/GRCm38.primary_assembly.genome.fa"
genome_region_bed = "/data/database/GRCm38/GENCODE/Genes_ensembl.bed"
#
promoter_bed = "/data5/galaxy/project/data/promoter/mouse/2k_100/genes_promoters.bed"
#
result_dir = "/data5/galaxy/project/DNMT1_KO/CpG_analysis"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = os.path.join(result_dir, "up_peak_genes_vs_others.txt")
###############################################################################################


def main():
    peak_genes = get_up_peak_genes()
    total_genes = filter_gene_expression(peak_genes)
    df_total_bed, df_peak_bed = get_bed(total_genes), get_bed(peak_genes)
    str_total_pro_bed, str_peak_pro_bed = get_promoter_bed(df_total_bed), get_promoter_bed(df_peak_bed)
    total_scores, peak_scores = enrich_CpG_in_promoter(str_total_pro_bed), enrich_CpG_in_promoter(str_peak_pro_bed)
    write_data_into_file(peak_scores, total_scores)
    result = stats.mannwhitneyu(peak_scores, total_scores)
    print(np.median(peak_scores), np.median(total_scores))
    print(result)


def filter_gene_expression(peak_genes):
    df = pd.read_csv(expression_file, index_col=0)
    df.index = [x.split(".")[0] for x in df.index]
    df_DTK_In = df.loc[:, "DTK_In"]
    df_filter = pd.DataFrame(df_DTK_In[df_DTK_In > 0])
    total_genes = df_filter.index
    #
    total_genes = [gene for gene in total_genes if not gene in peak_genes]
    return total_genes


def get_up_peak_genes():
    peak_genes = []
    with open(up_peak_gene, 'r') as f:
        for line in f.readlines():
            peak_genes.append(line.strip().split(".")[0])
    return peak_genes


def get_bed(genes):
    df = pd.DataFrame({"name": genes})
    df_bed = pd.read_table(genome_region_bed, sep="\t", header=None, names=["chr", "a", "b", "name", "c", "d"])
    df_overlap = pd.merge(df, df_bed, how="left")
    df_bed = df_overlap[["chr", "a", "b", "name", "c", "d"]].dropna(how="any").drop_duplicates()
    df_bed["a"], df_bed["b"] = df_bed["a"].astype(int), df_bed["b"].astype(int)
    df_gene_bed = df_bed.sort_values(["chr", "a"])
    return df_gene_bed


def get_promoter_bed(df_gene_bed):
    df_gene_bed.columns = ["a", "b", "c", "name", "e", "f"]
    df_promoter = pd.read_table(promoter_bed, sep="\t", header=None,
                                names=["chr", "start", "end", "name", "score", "strand"])
    df_merge = pd.merge(df_gene_bed, df_promoter, on="name", how="left")
    for col in ["a", "b", "c", "e", "f"]:
        del df_merge[col]
    df_bed = df_merge.dropna(how="any").drop_duplicates()
    df_bed["start"] = df_bed["start"].astype(int)
    df_bed["end"] = df_bed["end"].astype(int)
    df_bed["score"] = df_bed["score"].astype(int)
    df_promoter_bed = df_bed[["chr", "start", "end", "name", "score", "strand"]]
    str_promoter = str(df_promoter_bed.to_string(header=False, index=False, col_space=4)).strip()
    str_promoter = "\n".join(["\t".join(x.split()) for x in str_promoter.split("\n")])
    return str_promoter


def enrich_CpG_in_promoter(str_promoter):
    str_seq_list, score_list = get_sequence_from_bed(str_promoter).split("\n"), []
    for i in range(0, len(str_seq_list)-1, 2):
        title, seq = str_seq_list[i], str_seq_list[i+1]
        score_list.append(calculate_CpG_density(seq))
    return score_list


def get_sequence_from_bed(str_promoter):
    command = "bedtools getfasta -fi %s -bed stdin" % reference_genome
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    str_sequence = sub_p.communicate(str_promoter)[0]
    return str_sequence


def calculate_CpG_density(sequence):
    seq = sequence.lower()
    cg_num, c_num, g_num = len(re.findall("cg", seq)), len(re.findall("c", seq)), len(re.findall("g", seq))
    score = cg_num / (math.pow((c_num + g_num) / 2.0, 2) / (len(seq) * 1.0))
    return score


def write_data_into_file(pos_scores, neg_scores):
    with open(result_file, 'w') as fw:
        for score in pos_scores:
            fw.write("up peak genes\t%f\n" % score)
        for score in neg_scores:
            fw.write("other genes(in DTK In)\t%f\n" % score)


if __name__ == '__main__':
    main()