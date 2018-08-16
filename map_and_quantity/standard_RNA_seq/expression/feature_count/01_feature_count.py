#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
# @ function: only aimed at strand_specific data, now!
# when running featureCounts in its default setting, multi-mapping reads will be excluded.
# The '-M' option is only useful if you want to count those multi-mapping reads
import os
import glob
from multiprocessing import Pool


#########################################################################################
species, gff, bam_dir, result_dir = "mouse", "", "", ""
if species == "mouse":
    gff = "/data/database/GRCm38/GENCODE/gencode.vM15.annotation.gff3"
    bam_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/unique_bam"
    result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/expression/featureCount"
elif species == "human":
    gff = "/data/database/GRCh38/GENCODE/gencode.v27.annotation.gff3"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
result_file = "%s/gene_counts.txt" % result_dir
thread_num = 30
################################################################


def main(library_type):
    # base_dir = "/data5/galaxy/project/methyl_m6a/data/diff_m6a_peak/macs2"
    # tissue_list = [x.split("/")[-2] for x in glob.glob("%s/*/*.xls" % base_dir)]
    # pair_list_list = list(combinations(tissue_list, 2))
    # pool = Pool(processes=2)
    bam_list = [x for x in glob.glob("%s%s*_sorted_unique.bam" % (bam_dir, os.sep)) if "_IP_" not in x]
    bam_files = " ".join(bam_list)
    feature_count(library_type, bam_files)
    # for pair_list in pair_list_list:
    #     single_run(library_type, pair_list)
    # pool.close()
    # pool.join()


# def single_run(library_type, ):
    # tissue_1, tissue_2 = pair_list[0], pair_list[1]
    # print(tissue_1, tissue_2)
    # bam_dir = "/data5/galaxy/project/unique_bam/input/rename_bam/"
    # tissue_1_bams, tissue_2_bams = glob.glob("%s/%s*.bam" % (bam_dir, tissue_1)), glob.glob("%s/%s*.bam" % (bam_dir, tissue_2))
    # bam_files = " ".join(tissue_1_bams + tissue_2_bams)
    # print(bam_files)
    # result_dir = "/data5/galaxy/project/methyl_m6a/data/featureCount/%s_vs_%s" % (tissue_1, tissue_2)
    # parameter_setting(library_type, bam_files, result_dir)


def feature_count(library_type, bam_files):
    if library_type == "strand-specific":
        os.system("featureCounts -T %d -s 2 -p -a %s -o %s %s" % (thread_num, gff, result_file, bam_files))
    elif library_type == "normal":
        os.system("featureCounts -T %d -s 0 -a %s -o %s %s" % (thread_num, gff, result_file, bam_files))
    else:
        print("Can't configure library type!")


if __name__ == "__main__":
    # for lib_type in ["strand-specific"]:     # "normal, strand-specific"
    #     parameter_setting(lib_type)
    main("strand-specific")
