#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
# function: generate empirical null distribution; random 100 control region;
# one tissue, one distribution


import os
import glob
import pick_size_chr_match
from multiprocessing import Pool


#######################################################################################
ip_dir = "/data5/galaxy/project/data/ts_ip"
input_dir = "/data5/galaxy/project/bam2bed_input/merged_data/consecutive_bed"
cycle_100_dir = "/data5/galaxy/project/tf_analysis/bg_distribution/100cycle_input"
reference_genome = "/data4/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
motif_file = "/data4/database/tf_motif/motif_matrix.txt"
result_dir = "/data5/galaxy/project/tf_analysis/bg_distribution"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
count_result_file = os.path.join(result_dir, "background_distribution.txt")
#########################################################################################


def merge_all_peaks(bed_list, result_bed):
    bed_files = " ".join(bed_list)
    os.system("cat %s | sort -k1,1 -k2,2n - | bedtools merge -i - > %s" % (bed_files, result_bed))


def pro_process():
    input_list = glob.glob("%s/*.bed" % input_dir)
    print(input_list)
    ip_list = glob.glob("%s/*.bed" % ip_dir)
    print(ip_list)
    input_bed = os.path.join(result_dir, "input.bed")
    ip_bed = os.path.join(result_dir, "ip.bed")
    merge_all_peaks(input_list, input_bed)
    print("merge all input beds!")
    merge_all_peaks(ip_list, ip_bed)
    print("merge all ip beds!")
    return input_bed, ip_bed


def drawn_100_control(input_bed, ip_bed):
    pick_size_chr_match.main_run(input_bed, ip_bed, cycle_100_dir)


def transform_bed_to_fa(bed_file):
    fasta_file = "%s.fa" % bed_file.split(".bed")[0]
    os.system("bedtools getfasta -fi %s -bed %s -fo %s" % (reference_genome, bed_file, fasta_file))
    return fasta_file


def fasta_get_markov(fasta_file):
    markov_file = os.path.join(result_dir, "%s.txt" % os.path.basename(fasta_file).split(".")[0])
    os.system("fasta-get-markov %s %s" % (fasta_file, markov_file))
    return markov_file


def fimo_transcript_factor_num(bed_file, out_dir):
    fasta_file = transform_bed_to_fa(bed_file)
    markov_file = fasta_get_markov(fasta_file)
    os.system("fimo --bgfile %s --oc %s %s %s" % (markov_file, out_dir, motif_file, fasta_file))


def calculate_enrich_ratio():
    result_list = []
    fimo_list = glob.glob("%s/*/fimo.txt" % result_dir)
    ip_list = [i_result for i_result in fimo_list if i_result.split("/")[-2] == "ip"]
    input_list = [i_result for i_result in fimo_list if i_result.split("/")[-2] != "ip"]
    ip_count, input_count_list = statistic_motif_num(ip_list[0]), []
    for i_input in input_list:
        input_count_list.append(statistic_motif_num(i_input))
    for count in input_count_list:
        result_list.append(float(ip_count) / float(count))
    with open(count_result_file, 'w') as fw:
        for line in result_list:
            fw.write("%s\n" % line)


def statistic_motif_num(in_file):
    with open(in_file, 'r') as f:
        count = len(f.readlines()) - 1
    return count


def search_tf_motif():
    """
    input_bed, ip_bed = pro_process()
    drawn_100_control(input_bed, ip_bed)
    input_list = glob.glob("%s/*.bed" % cycle_100_dir)
    total_list = input_list + [ip_bed]
    pool = Pool(processes=20)
    for i_bed in total_list:
        print(i_bed)
        out_dir = os.path.join(result_dir, os.path.basename(i_bed).split(".bed")[0])
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        pool.apply_async(fimo_transcript_factor_num, (i_bed, out_dir))
        # fimo_transcript_factor_num(i_bed, out_dir)
    pool.close()
    pool.join()
    """
    calculate_enrich_ratio()


if __name__ == "__main__":
    search_tf_motif()
