#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool

fq_dir, ribo_gtf, ribo_index, ref_index, result_dir, bowtie2_result_dir, hisat2_result_dir = "", "", "", "", "", "", ""
species = "mouse"
if species == "mouse":
    fq_dir = "/data4/galaxy/project/plot_what"
    # ribo_genome = "/data4/database/mm10/mm10_Repeats_rRNA.fa"
    # ribo_gtf = "/data4/database/mm10/mm10_Repeats_rRNA.gtf"
    ribo_index = "/data4/database/mm10/rRNA"
    ref_index = "/data4/database/mm10/mm10"
    splice_sites = "/data4/database/mm10/splicesites.txt"
    bowtie2_result_dir = "/data4/galaxy/project/plot_what/map_result/01bowtie2_result"
    hisat2_result_dir = "/data4/galaxy/project/plot_what/map_result/02hisat2_result"
    # result_dir = "/data4/galaxy/project/plot_what/hisat_result"
elif species == "human":
    fq_dir = "/data3/xs/original_data/d0905-P101SC17080761-01-B1-RRL-N17/data_release/clean_data"
    ribo_index = "/data2/galaxy/project/correlation/rm_ribo/hg38_ribo_database/hg38_ribo_bowtie2_index/hg38"
    ref_index = "/data2/galaxy/project/correlation/rm_ribo/hisat2_index/hg38"
    splice_sites = "/data4/database/hg38/splicesites.txt"
    bowtie2_result_dir = "/data4/galaxy/project/new_analysis_m6a/map_result/01bowtie2_result"
    hisat2_result_dir = "/data4/galaxy/project/new_analysis_m6a/map_result/02hisat2_result"
for i_result_dir in [bowtie2_result_dir, hisat2_result_dir]:
    if not os.path.exists(i_result_dir):
        os.makedirs(i_result_dir)
each_thread_num = 2
#####################################


def get_all_fasta():
    fq_files = glob.glob("%s/*_1.clean.fq.gz" % fq_dir)
    return fq_files


def bowtie2_ribosome(fq_1, fq_2):
    prefix = os.path.basename(fq_1).split("_")[0]
    mapped_fq = os.path.join(bowtie2_result_dir, "%s_mapped.fq" % prefix)
    unmapped_fq = os.path.join(bowtie2_result_dir, "%s_unmapped.fq" % prefix)
    ribo_sam = os.path.join(bowtie2_result_dir, "%s_ribo.sam" % prefix)
    os.system("bowtie2 --fr --al-conc %s --un-conc %s -p %i -x %s -1 %s -2 %s -S %s" % (mapped_fq, unmapped_fq, each_thread_num, ribo_index, fq_1, fq_2, ribo_sam))
    unmap_1, unmap_2 = unmapped_fq.replace("_unmapped.fq", "_unmapped.1.fq"), unmapped_fq.replace("_unmapped.fq", "_unmapped.2.fq")
    return unmap_1, unmap_2


def hisat2_each(fq_1, fq_2):
    prefix = os.path.basename(fq_1).split("_")[0]
    out_sam = os.path.join(hisat2_result_dir, "%s.sam" % prefix)
    unmapped_fq = os.path.join(hisat2_result_dir, "%s_unmapped.fq" % prefix)
    log_file = os.path.join(hisat2_result_dir, "%s.log" % prefix)
    os.system("hisat2 --known-splicesite-infile %s -x %s -1 %s -2 %s --dta -t -p %d --rna-strandness RF --un-conc %s -S %s 2> %s" % (splice_sites, ref_index, fq_1, fq_2, each_thread_num, unmapped_fq, out_sam, log_file))


def two_step_map(fq_1, fq_2):
    unmap_1, unmap_2 = bowtie2_ribosome(fq_1, fq_2)
    hisat2_each(unmap_1, unmap_2)


def map_run():
    fq_1_files = get_all_fasta()
    pool = Pool(processes=len(fq_1_files))
    for fq_1 in fq_1_files:
        fq_2 = fq_1.replace("_1.clean.fq.gz", "_2.clean.fq.gz")
        pool.apply_async(two_step_map, (fq_1, fq_2))
    pool.close()
    pool.join()

# fw.write("samtools view -b -q 5 -o %s %s\n" % (rmdup_bam, unique_bam))
# eRNA long [50, 2000]


if __name__ == "__main__":
    map_run()
