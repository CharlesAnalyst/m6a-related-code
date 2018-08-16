#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os

#######################################################################
species = "human"
each_thread_num = 35
ref_index = "/data4/database/GRCh38/GENCODE/GRCh38"
splice_sites = "/data4/database/GRCh38/GENCODE/splice_sites.txt"
ribo_index = "/data4/database/hg38_ribo_database/hg38_ribo_bowtie2_index/hg38"
bowtie2_result_dir = "/data5/galaxy/project/map_and_quantity/input/01bowtie2_result"
hisat2_result_dir = "/data5/galaxy/project/map_and_quantity/input/02hisat2_result"
postfix = "_1.clean.fq.gz"
#########
picard = "java -jar /home/xiaoshan/bin/picard.jar"
tmp_dir = "/data5/galaxy/project"
contain_word = ".sam"
result_dir = "/data5/galaxy/project/unique_bam/input"
#####################################


def get_all_fq_1():
    fq_1_files = ["/data3/xs/original_data/d0905-P101SC17080761-01-B1-RRL-N17/data_release/clean_data/brain_combine/s113-S1-brain_H7YM_LX_1.clean.fq.gz"]
    return fq_1_files


def get_prefix(fq):
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


def hisat2_pair_end(fq_1, fq_2):
    prefix = os.path.basename(fq_1).split("_unmapped.1.fq")[0]
    out_sam = os.path.join(hisat2_result_dir, "%s.sam" % prefix)
    unmapped_fq = os.path.join(hisat2_result_dir, "%s_unmapped.fq" % prefix)
    log_file = os.path.join(hisat2_result_dir, "%s.log" % prefix)
    os.system("hisat2 --rna-strandness RF --known-splicesite-infile %s -x %s -1 %s -2 %s --dta -t -p %d --un-conc %s -S %s 2> %s" % (splice_sites, ref_index, fq_1, fq_2, each_thread_num, unmapped_fq, out_sam, log_file))
    return out_sam


def get_prefix_from_sam(i_sam):
    prefix = os.path.basename(i_sam).split(contain_word)[0]
    return prefix


def get_prefix_from_sorted_bam(i_bam):
    prefix = os.path.basename(i_bam).split(".bam")[0]
    return prefix


def sam2bam_sort(i_sam):
    prefix = get_prefix_from_sam(i_sam)
    unsorted_bam = "%s/%s_unsorted.bam" % (result_dir, prefix)
    os.system("samtools view -@ %d -bS %s > %s" % (each_thread_num, i_sam, unsorted_bam))
    sorted_bam = "%s/%s_sorted" % (result_dir, prefix)
    os.system("samtools sort -@ %d %s %s" % (each_thread_num, unsorted_bam, sorted_bam))
    os.system("rm %s" % unsorted_bam)
    return sorted_bam


def rmdup_unique_multi_merge(i_bam):
    prefix = get_prefix_from_sorted_bam(i_bam)
    rmdup_bam = os.path.join(result_dir, "%s_rmdup.bam" % prefix)
    metric_file = os.path.join(result_dir, "%s_out.metrics" % prefix)
    os.system("%s MarkDuplicates I=%s O=%s REMOVE_DUPLICATES=true CREATE_INDEX=true TMP_DIR=%s VALIDATION_STRINGENCY=SILENT M=%s" % (picard, i_bam, rmdup_bam, tmp_dir, metric_file))
    unique_bam = os.path.join(result_dir, "%s_unique.bam" % prefix)
    os.system("samtools view -@ %d -b -q 30 -o %s %s" % (each_thread_num, unique_bam, rmdup_bam))
    os.system("samtools index %s" % unique_bam)
    # os.system("rm %s" % rmdup_bam)


def main_run():
    fq_1 = get_all_fq_1()[0]
    fq_2 = fq_1.replace("_1.", "_2.")
    unmap_1, unmap_2 = bowtie2_pair_end(fq_1, fq_2)
    out_sam = hisat2_pair_end(unmap_1, unmap_2)
    sorted_bam_prefix = sam2bam_sort(out_sam)
    sorted_bam = sorted_bam_prefix + ".bam"
    rmdup_unique_multi_merge(sorted_bam)


if __name__ == "__main__":
    main_run()
