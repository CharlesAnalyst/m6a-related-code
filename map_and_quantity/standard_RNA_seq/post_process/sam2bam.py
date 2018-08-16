#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool

################################################################
picard = "java -jar /home/xiaoshan/bin/picard.jar"
sam_dir = ""
# bam_type = "SE"  # bam_type: SE or PE
sam_list = glob.glob("%s/*.sam" % sam_dir)
result_dir = ""
contain_word = ""
total_thread_num = 30
thread_num = 15
tmp_dir = "/data5/galaxy/project"
################################################################


def parameter_setting(data_type):
    global sam_dir, contain_word, result_dir, sam_list, thread_num
    if data_type == "ip":
        sam_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/map_result/02hisat2_result"
        contain_word = ".sam"
        sam_list = glob.glob("%s/*%s" % (sam_dir, contain_word))
        # thread_num = total_thread_num / len(sam_list)
        result_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/uniq_bam"
    elif data_type == "input":
        sam_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/02hisat2_result"
        contain_word = ".sam"
        sam_list = glob.glob("%s/*%s" % (sam_dir, contain_word))
        thread_num = total_thread_num / len(sam_list)
        result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/unique_bam"
    else:
        print("Can't recognise the data type!")
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)


def get_prefix_from_sam(i_sam):
    prefix = os.path.basename(i_sam).split(contain_word)[0]
    return prefix


def get_prefix_from_sorted_bam(i_bam):
    prefix = os.path.basename(i_bam).split(".bam")[0]
    return prefix


def sam2bam_sort(i_sam):    # samtools sort -@ 8 -o HBR_Rep3.bam HBR_Rep3.sam
    prefix = get_prefix_from_sam(i_sam)
    unsorted_bam = "%s/%s_unsorted.bam" % (result_dir, prefix)
    os.system("samtools view -@ %d -bS %s > %s" % (thread_num, i_sam, unsorted_bam))
    #"samtools view -@ 30 -bS /data5/galaxy/project/map_and_quantity/input/02hisat2_result/s113-brain_H7YMKCCXY_LX.sam > s113-brain_H7YMKCCXY_LX_unsorted.bam"
    sorted_bam = "%s/%s_sorted" % (result_dir, prefix)
    os.system("samtools sort -@ %d %s %s" % (thread_num, unsorted_bam, sorted_bam))
    os.system("rm %s" % unsorted_bam)
    return sorted_bam


def rmdup_unique_multi_merge(i_bam):
    prefix = get_prefix_from_sorted_bam(i_bam)
    rmdup_bam = os.path.join(result_dir, "%s_rmdup.bam" % prefix)
    metric_file = os.path.join(result_dir, "%s_out.metrics" % prefix)
    os.system("%s MarkDuplicates I=%s O=%s REMOVE_DUPLICATES=true CREATE_INDEX=true TMP_DIR=%s VALIDATION_STRINGENCY=SILENT M=%s" % (picard, i_bam, rmdup_bam, tmp_dir, metric_file))
    # pick out unique mapped reads
    unique_bam = os.path.join(result_dir, "%s_unique.bam" % prefix)
    os.system("samtools view -@ %d -b -q 20 -o %s %s" % (thread_num, unique_bam, rmdup_bam))
    # pick out multi mapped reads
    # multi_bam = os.path.join(result_dir, "%s_multi.bam" % prefix)
    # os.system("fgrep -vw NH:i:1 %s > %s" % (rmdup_bam, multi_bam))
    # merge unique + multi bam
    # merge_bam = os.path.join(result_dir, "%s_merged.bam" % prefix)
    # os.system("samtools merge -nr -@ %d %s %s %s" % (thread_num, merge_bam, unique_bam, multi_bam))
    # build bam index
    # for i_bam in [unique_bam, multi_bam, merge_bam]:
    #    os.system("samtools index %s" % i_bam)
    #
    os.system("samtools index %s" % unique_bam)
    if os.path.exists(unique_bam):
        os.system("rm %s" % rmdup_bam)


def two_step_process(i_sam):
    sorted_bam_prefix = sam2bam_sort(i_sam)
    sorted_bam = sorted_bam_prefix + ".bam"
    rmdup_unique_multi_merge(sorted_bam)


if __name__ == "__main__":
    # for i_type in ["ip", "input"]:
    for i_type in ["input"]:
        parameter_setting(i_type)
        for sam in sam_list:
            two_step_process(sam)
        """"
        pool = Pool()
        for sam in sam_list:
            pool.apply_async(two_step_process, (sam,))
        pool.close()
        pool.join()
"""