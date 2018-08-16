#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
# @ function: only aimed at strand_specific data, now!
# when running featureCounts in its default setting, multi-mapping reads will be excluded.
# The '-M' option is only useful if you want to count those multi-mapping reads
import os
import glob

###############################################################
gff = "/data/database/GRCm38/GENCODE/gencode.vM15.annotation.gff3"
thread_num = 44
# strand_specific = "yes"  # can't alter
################################################################


def parameter_setting(data_type):
    if data_type == "PE":
        bam_dir = "/data5/galaxy/project/unique_bam/input"
        bam_files = " ".join(glob.glob("%s/*_sorted_unique.bam" % bam_dir))
        result_dir = "/data5/galaxy/project/expression/featureCounts/input"
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        result_file = "%s/gene_counts.txt" % result_dir
        os.system("featureCounts -T %d -s 2 -p -a %s -o %s %s" % (thread_num, gff, result_file, bam_files))
    elif data_type == "SE":
        bam_dir = "/data5/galaxy/project/DNMT1_KO/unique_bam"
        bam_files = " ".join(glob.glob("%s/*_sorted_unique.bam" % bam_dir))
        result_dir = "/data5/galaxy/project/DNMT1_KO/featureCounts"
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        result_file = "%s/gene_counts.txt" % result_dir
        os.system("featureCounts -T %d -s 0 -a %s -o %s %s" % (thread_num, gff, result_file, bam_files))
    else:
        print("Can't configure data type!")


if __name__ == "__main__":
    for i_type in ["SE"]:
        parameter_setting(i_type)