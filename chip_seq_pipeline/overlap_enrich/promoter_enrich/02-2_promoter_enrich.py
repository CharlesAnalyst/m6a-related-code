#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import subprocess
import pandas as pd
from scipy import stats
from multiprocessing import Pool

###############################################################################
gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/gene_bed"
gene_list = glob.glob("%s/*.bed" % gene_dir)
gene_bed_dict = {}
for g_bed in gene_list:
    tis_name = os.path.basename(g_bed).split(".bed")[0].lower()
    gene_bed_dict[tis_name] = g_bed
#########
m6a_bed_dir = "/data5/galaxy/project/data/total_m6a_peak"
m6a_list = glob.glob("%s/*.bed" % m6a_bed_dir)
m6a_gene_dir = "/data5/galaxy/project/promoter_TF_enrich/m6a_gene"
#
# tfbs_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_hg19_to_GRCh38"
stat1_dir_one = "/data3/xs/tissue_m6a/2018.1/TF/ENCODE_PEAKS/STAT1"
stat1_dir_two = "/data3/xs/tissue_m6a/2018.1/TF/ENCODE_PEAKS/STAT1_hela"
thap11_dir = "/data5/galaxy/project/tf_analysis/TFBS_analysis/data/THAP11"
tfbs_list = []
for data_dir in [stat1_dir_one, stat1_dir_two, thap11_dir]:
    tfbs_list += glob.glob("%s/*.bed" % data_dir)
print(len(tfbs_list))
#
result_dir = "/data5/galaxy/project/promoter_TF_enrich/part_TF/tissue_enr_result"
gene_promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/gene_promoter_bed"
promoter_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/promoter_bed"
for i_dir in [result_dir, gene_promoter_dir, promoter_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)
promoter_bed, gene_promoter_bed, gene_dict, promoter_dict = "", "", {}, {}
####################################################################################################


def get_name_bed_dict():
    result_dict = {}
    with open(promoter_bed, 'r') as f:
        for line in f.readlines():
            name = line.split("\t")[3]
            result_dict[name] = line.strip()
    return result_dict


def get_bed_name_dict():
    result_dict = {}
    with open(gene_promoter_bed, 'r') as f:
        for line in f.readlines():
            name = line.split("\t")[3]
            result_dict[line.strip()] = name
    return result_dict


# TSS upstream 2kb
def get_gene_promoter(gene_bed):
    gene_promoter_bed = os.path.join(gene_promoter_dir, os.path.basename(gene_bed))
    df_gene = pd.read_table(gene_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "a", "strand"])
    df_gene_positive = df_gene[df_gene.loc[:, "strand"] == "+"]
    df_gene_negative = df_gene[df_gene.loc[:, "strand"] == "-"]
    #
    df_gene_negative.loc[:, "promoter_end"] = df_gene_negative["end"] + 2000
    df_gene_negative_bed = df_gene_negative[["chr", "start", "promoter_end", "name", "a", "strand"]]
    df_gene_negative_bed.columns = ["chr", "start", "end", "name", "a", "strand"]
    #
    df_gene_positive.loc[:, "promoter_start"] = df_gene_positive["start"] - 2000
    df_gene_positive_bed = df_gene_positive[["chr", "promoter_start", "end", "name", "a", "strand"]]
    df_gene_positive_bed.columns = ["chr", "start", "end", "name", "a", "strand"]
    #
    df_promoter_list = [df_gene_positive_bed, df_gene_negative_bed]
    df_promoter = pd.concat(df_promoter_list)
    df_promoter_sort = df_promoter.sort_values(["chr", "start"])
    #
    df_promoter_sort = df_promoter_sort[df_promoter_sort["start"] > 0]
    df_promoter_sort.to_csv(gene_promoter_bed, sep="\t", header=None, index=False)
    return gene_promoter_bed


# TSS upstream 2kb
def get_promoter(gene_bed):
    promoter_bed = os.path.join(promoter_dir, os.path.basename(gene_bed))
    df_gene = pd.read_table(gene_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "a", "strand"])
    df_gene_positive = df_gene.loc[df_gene.loc[:, "strand"] == "+", :]
    df_gene_negative = df_gene.loc[df_gene.loc[:, "strand"] == "-", :]
    df_gene_negative["promoter_end"] = df_gene_negative.loc[:, "end"] + 2000
    # sys.exit(1)
    df_gene_negative.loc[:, "promoter_end"] = df_gene_negative.loc[:, "end"] + 2000
    # print("bbbb")
    df_gene_negative_bed = df_gene_negative.loc[:, ["chr", "end", "promoter_end", "name", "a", "strand"]]
    # print(df_gene_negative_bed.head())
    df_gene_negative_bed.columns = ["chr", "start", "end", "name", "a", "strand"]
    #
    df_gene_positive.loc[:, "promoter_start"] = df_gene_positive.loc[:, "start"] - 2000
    df_gene_positive_bed = df_gene_positive[["chr", "promoter_start", "start", "name", "a", "strand"]]
    # print(df_gene_positive_bed.head())
    df_gene_positive_bed.columns = ["chr", "start", "end", "name", "a", "strand"]
    #
    df_promoter_list = [df_gene_positive_bed, df_gene_negative_bed]
    df_promoter = pd.concat(df_promoter_list)
    df_promoter_sort = df_promoter.sort_values(["chr", "start"])
    df_promoter_sort = df_promoter_sort[df_promoter_sort["start"] > 0]
    df_promoter_sort.to_csv(promoter_bed, sep="\t", header=None, index=False)
    # print(df_promoter_sort.head())
    # """
    return promoter_bed


def get_promoter_according_gene(gene_bed_list):
    bed_list = [promoter_dict[gene_dict[gene]] for gene in gene_bed_list if gene != ""]
    return bed_list


def class_gene_according_m6a(m6a_bed):
    with open(gene_promoter_bed, 'r') as f:
        total_genes = [x.strip() for x in f.readlines()]
    command = "bedtools intersect -a %s -b %s -wa | uniq" % (gene_promoter_bed, m6a_bed)
    m6a_genes = str(subprocess.check_output(command, shell=True)).split("\n")
    non_m6a_genes = [x for x in total_genes if x not in m6a_genes]
    # print("m6a genes && non_m6a_genes")
    # print(len(m6a_genes), len(non_m6a_genes))
    return m6a_genes, non_m6a_genes


# only promoter
def class_gene_according_chip(in_txt, tfbs_bed):
    command = "bedtools intersect -a stdin -b %s -F 0.5 -wa | uniq" % tfbs_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    tf_genes_num = len(str(sub_p.communicate(in_txt)[0]).split("\n"))
    non_tf_gene_num = len(in_txt.split("\n")) - tf_genes_num
    return tf_genes_num, non_tf_gene_num


def run_fisher_test(a, b, c, d):
    print("################################################")
    print("# ###              m6a_gene      non_m6a_gene")
    print("# ### tf_gene         %s               %s" % (a, b))
    print("# ### non_tf_gene     %s               %s" % (c, d))
    print("################################################")
    oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
    print(a+b+c+d)
    return oddsratio, pvalue


# promoter + gene body
def class_gene_to_table(m6a_bed, tfbs_bed):
    tissue = os.path.basename(m6a_bed).split(".bed")[0].lower()
    tfbs_name = os.path.basename(tfbs_bed).split(".bed")[0].upper()
    #
    m6a_genes, non_m6a_genes = class_gene_according_m6a(m6a_bed)
    # transform gene to promoter
    m6a_genes = get_promoter_according_gene(m6a_genes)
    non_m6a_genes = get_promoter_according_gene(non_m6a_genes)
    #
    m6a_tf_num, m6a_non_tf_num = class_gene_according_chip("\n".join(m6a_genes), tfbs_bed)
    non_m6a_tf_num, non_m6a_non_tf_num = class_gene_according_chip("\n".join(non_m6a_genes), tfbs_bed)
    #
    oddsratio, pvalue = run_fisher_test(m6a_tf_num, non_m6a_tf_num, m6a_non_tf_num, non_m6a_non_tf_num)
    print("%s\t%s\t%f\t%f" % (tissue, tfbs_name, oddsratio, pvalue))
    return str("%s\t%s\t%f\t%f\n" % (tissue, tfbs_name, oddsratio, pvalue))


def write2file(in_list, result_file):
    with open(result_file, 'w') as fw:
        fw.writelines(in_list)


def process_by_each_tissue(tissue_bed, gene_bed):
    results = []
    global promoter_bed, gene_promoter_bed, gene_dict, promoter_dict
    promoter_bed, gene_promoter_bed = get_promoter(gene_bed), get_gene_promoter(gene_bed)
    """
    gene_dict, promoter_dict = get_bed_name_dict(), get_name_bed_dict()
    pool = Pool()
    for tfbs_bed in tfbs_list:
        result = pool.apply_async(class_gene_to_table, (tissue_bed, tfbs_bed))
        # result = class_gene_to_table(tissue_bed, tfbs_bed)
        results.append(result)
        # break
    pool.close()
    pool.join()
    tissue_name = os.path.basename(tissue_bed).split(".bed")[0]
    result_file = os.path.join(result_dir, "%s.txt" % tissue_name)
    result_list = [x.get() for x in results]
    write2file(result_list, result_file)
"""


if __name__ == '__main__':
    for i_bed in m6a_list:
        print(i_bed)
        g_bed = gene_bed_dict[os.path.basename(i_bed).split(".bed")[0].lower()]
        process_by_each_tissue(i_bed, g_bed)
