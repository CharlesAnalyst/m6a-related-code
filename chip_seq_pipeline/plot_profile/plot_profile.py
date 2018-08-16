#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
from multiprocessing import Pool

tf_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean"
stat1, thap11, p53, rela = "%s/STAT1.bed" % tf_dir, "%s/THAP11.bed" % tf_dir, "%s/P53.bed" % tf_dir, "%s/RELA.bed" % tf_dir
tf_list = [stat1, thap11, p53, rela]
bw_dir = "/data5/galaxy/project/mettl3_enrich/clean_bam/plot_profile/current_version"
ip_bw, input_bw = "%s/ip_coverage.bw" % bw_dir, "%s/input_coverage.bw" % bw_dir
bw_list = [ip_bw, input_bw]
result_dir = bw_dir
# "nohup bamCoverage -b stat1_ip_unique.bam -o stat1_ip.bw --normalizeUsingRPKM &
#  --effectiveGenomeSize 2701495761
# computeMatrix reference-point -S %s -R %s -a 3000 -b 3000 -o %s
#
# computeMatrix reference-point --referencePoint center -S ip_coverage.bw input_coverage.bw -R .promoter_3_3bed -a 3000 -b 3000 -o mettl3_tss.gz
# plotProfile -m mettl3_tss.gz -out mettl3_tss_3kb.pdf --perGroup --plotTitle ""  # make one image per BED file instead of per bigWig file


def compute_matrix(bw_file, bed, result):
    os.system("computeMatrix reference-point -S %s -R %s -a 3000 -b 3000 -o %s" % (bw_file, bed, result))
    # nohup computeMatrix reference-point -S stat1_*bw -R /data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed -a 3000 -b 3000 --referencePoint center -p max -o stat1.gz &


def multi_compute_matrix():
    for bw in bw_list:
        print(bw)
        bw_prefix = os.path.basename(bw).split("_")[0]
        pool = Pool()
        for tf in tf_list:
            print(tf)
            tf_prefix = os.path.basename(tf).split(".bed")[0]
            result_file = os.path.join(result_dir, "%s_%s.gz" % (bw_prefix, tf_prefix))
            pool.apply_async(compute_matrix, (bw, tf, result_file))
        pool.close()
        pool.join()

#
def plot_profile(matrix_file):
    picture_file = matrix_file.replace(".gz", ".png")
    title = os.path.basename(matrix_file).split(".gz")[0]
    os.system("plotProfile -m %s -out %s --plotTitle %s" % (matrix_file, picture_file, title))


def plot_heatmap(matrix_file):
    picture_file = matrix_file.replace(".gz", "_heatmap.png")
    # title = os.path.basename(matrix_file).split(".gz")[0]
    os.system("plotHeatmap -m %s -out %s " % (matrix_file, picture_file))


def multi_plot_profile():
    matrix_list = glob.glob("%s/*.gz" % bw_dir)
    pool = Pool()
    for matrix in matrix_list:
        # pool.apply_async(plot_profile, (matrix, ))
        pool.apply_async(plot_heatmap, (matrix,))
    pool.close()
    pool.join()


# nohup computeMatrix scale-regions -R /data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean/THAP11.bed /data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean/P53.bed /data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean/STAT1.bed /data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean/MYC.bed -S /data4/shijunfang_1/mettl3_chip_hela/clean_data/mettl3.bw -b 3000 -a 3000 --regionBodyLength 5000 --skipZeros -o mettl3_ip_total.gz --outFileNameMatrix mettl3_ip_total.tab --outFileSortedRegions mettl3_ip_total.bed &


if __name__ == '__main__':
    multi_plot_profile()
