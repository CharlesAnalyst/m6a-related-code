#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
from multiprocessing import Pool


anno_type = "3UTR"    # 3UTR
annotate_dir = "/data5/galaxy/project/data/total_m6a_peak/annotate_peak/%s" % anno_type
peak_list = glob.glob("%s/*.bed" % annotate_dir)

gene_dir = "/data5/galaxy/project/promoter_TF_enrich/data/total_gene/gene_bed"
gene_list = glob.glob("%s/*.bed" % gene_dir)


result_dir = "/data5/galaxy/project/promoter_TF_enrich/m6a_vs_free/3-5_UTR/%s/gene_bed" % anno_type
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def main():
    pool = Pool()
    for peak in peak_list:
        gene_bed = os.path.join(gene_dir, os.path.basename(peak))
        # select_genes(gene_bed, peak)
        pool.apply_async(select_genes, (gene_bed, peak))
    pool.close()
    pool.join()


def select_genes(gene_bed, UTR_bed):
    intersect_result_file = os.path.join(result_dir, os.path.basename(gene_bed))
    os.system("bedtools intersect -a {gene_bed} -b {UTR_bed} -wa > {intersect_result_file}"
              .format(gene_bed=gene_bed, UTR_bed=UTR_bed, intersect_result_file=intersect_result_file))


if __name__ == '__main__':
    main()