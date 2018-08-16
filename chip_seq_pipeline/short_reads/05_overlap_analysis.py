#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import permutation_test as p


up_peak_gene = "/data5/galaxy/project/DNMT1_KO/diff_peak/DTK_vs_CGR8_c3.0_cond1_anno.txt"
up_gene = "/data5/galaxy/project/DNMT1_KO/diff_expre/up_genes.txt"
result_dir = "/data5/galaxy/project/DNMT1_KO/overlap_analysis"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def main():
    up_peak_genes = read_file(up_peak_gene)
    up_genes = read_file(up_gene)
    plot_venn2(up_peak_genes, up_genes)
    pvalue = permutation_test(up_peak_genes, up_genes)
    print(pvalue)


def read_file(each_file):
    with open(each_file, 'r') as f:
        genes = f.readlines()
        gene_list = [x.strip().split(".")[0].upper() for x in genes]
    return gene_list


def plot_venn2(list_1, list_2):
    plt.figure(figsize=(4, 4))
    venn2([set(list_1), set(list_2)], set_labels=("up_peak_genes", "up_express_genes"))
    # plt.show()
    plt.savefig('%s/venn.pdf' % result_dir)


def permutation_test(list_1, list_2):
    p_value = p.permutation_test(list_1, list_2)
    return p_value


if __name__ == '__main__':
    main()