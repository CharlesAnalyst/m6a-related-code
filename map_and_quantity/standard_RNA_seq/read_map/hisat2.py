#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool

#######################################################################
species = "mouse"
rm_ribosome = "yes"
total_thread_num = 30
##########################
if species == "mouse":
    ref_index = "/data/database/GRCm38/GENCODE/GRCm38"
    if rm_ribosome == "yes":
        ribo_index = "/data/database/mm10/rRNA"
elif species == "human":
    ref_index = "/data/database/GRCh38/GENCODE/hisat2_index/grch38_snp_tran/genome_snp_tran"
    if rm_ribosome == "yes":
        ribo_index = "/data/database/hg38_ribo_database/hg38_ribo_bowtie2_index/hg38"
data_type, postfix, fq_dir, hisat2_result_dir, bowtie2_result_dir, each_thread_num = "", "", "", "", "", 1
#####################################


def parameter_setting(case_type):
    global data_type, postfix, fq_dir, hisat2_result_dir, bowtie2_result_dir
    if case_type == "ip":
        data_type = "SE"
        postfix = ".fastq"
        fq_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/m6a_level/clean_fq"
        bowtie2_result_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/map_result/01bowtie2_result"
        hisat2_result_dir = "/data5/galaxy/project/DNMT1_m6a/MEF_MES_compare/map_result/02hisat2_result"
    elif case_type == "input":
        data_type = "PE"
        postfix = "_1.clean.fq.gz"
        fq_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/fq_dir"
        bowtie2_result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/01bowtie2_result_2"
        hisat2_result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/02hisat2_result_2"
    for i_result_dir in [bowtie2_result_dir, hisat2_result_dir]:
        if not os.path.exists(i_result_dir):
            os.makedirs(i_result_dir)


def get_all_fq_1():
    fq_1_files = glob.glob("%s/*%s" % (fq_dir, postfix))
    global each_thread_num
    each_thread_num = total_thread_num / len(fq_1_files)
    # each_thread_num = total_thread_num
    return fq_1_files


def get_prefix(fq):  # get prefix from fq_1
    prefix = os.path.basename(fq).split(postfix)[0]
    return prefix


def bowtie2_pair_end(fq_1, fq_2):
    prefix = get_prefix(fq_1)
    mapped_fq = os.path.join(bowtie2_result_dir, "%s_mapped.fq" % prefix)
    unmapped_fq = os.path.join(bowtie2_result_dir, "%s_unmapped.fq" % prefix)
    ribo_sam = os.path.join(bowtie2_result_dir, "%s_ribo.sam" % prefix)
    os.system("bowtie2 --fr --al-conc %s --un-conc %s -p %i -x %s -1 %s -2 %s -S %s" % (mapped_fq, unmapped_fq, each_thread_num, ribo_index, fq_1, fq_2, ribo_sam))
    unmap_1, unmap_2 = unmapped_fq.replace("_unmapped.fq", "_unmapped.1.fq"), unmapped_fq.replace("_unmapped.fq", "_unmapped.2.fq")
    return unmap_1, unmap_2


def bowtie2_single_end(fq):
    prefix = get_prefix(fq)
    mapped_fq = os.path.join(bowtie2_result_dir, "%s_mapped.fq" % prefix)
    unmapped_fq = os.path.join(bowtie2_result_dir, "%s_unmapped.fq" % prefix)
    ribo_sam = os.path.join(bowtie2_result_dir, "%s_ribo.sam" % prefix)
    os.system("bowtie2 --al %s --un %s -p %i -x %s -U %s -S %s" % (mapped_fq, unmapped_fq, each_thread_num, ribo_index, fq, ribo_sam))
    return unmapped_fq


def hisat2_pair_end(fq_1, fq_2):
    prefix = os.path.basename(fq_1).split("_unmapped.1.fq")[0]
    out_sam = os.path.join(hisat2_result_dir, "%s.sam" % prefix)
    unmapped_fq = os.path.join(hisat2_result_dir, "%s_unmapped.fq" % prefix)
    log_file = os.path.join(hisat2_result_dir, "%s.log" % prefix)
    os.system("hisat2 --rna-strandness RF -x %s -1 %s -2 %s --dta -t -p %d --un-conc %s -S %s 2> %s" % (ref_index, fq_1, fq_2, each_thread_num, unmapped_fq, out_sam, log_file))


def hisat2_single_end(fq):
    prefix = os.path.basename(fq).split("_unmapped.fq")[0]
    out_sam = os.path.join(hisat2_result_dir, "%s.sam" % prefix)
    # unmapped_fq = os.path.join(hisat2_result_dir, "%s_unmapped.fq" % prefix)
    log_file = os.path.join(hisat2_result_dir, "%s.log" % prefix)
    if species == "mouse":
        os.system("hisat2 -x %s -U %s --dta -t -p %d -S %s 2> %s"
                  % (ref_index, fq, each_thread_num, out_sam, log_file))
    else:
        os.system("hisat2 -x %s -U %s --dta -t -p %d -S %s 2> %s" % (ref_index, fq, each_thread_num, out_sam, log_file))


def two_step_map(fq_1):
    if data_type == "PE":
        fq_2 = fq_1.replace("_1.", "_2.")
        unmap_1, unmap_2 = bowtie2_pair_end(fq_1, fq_2)
        hisat2_pair_end(unmap_1, unmap_2)
    else:
        unmapped_fq = bowtie2_single_end(fq_1)
        hisat2_single_end(unmapped_fq)


def map_run():
    fq_1_files = get_all_fq_1()
    pool = Pool(processes=len(fq_1_files))
    for fq_1 in fq_1_files:
        pool.apply_async(two_step_map, (fq_1, ))
        # two_step_map(fq_1)
    pool.close()
    pool.join()


if __name__ == "__main__":
    for i_type in ["input"]:
        parameter_setting(i_type)
        map_run()
