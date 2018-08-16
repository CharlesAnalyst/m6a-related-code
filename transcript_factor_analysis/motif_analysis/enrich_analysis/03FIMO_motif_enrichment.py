#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import datetime
from multiprocessing import Pool


#######################################################################################
ip_dir = "/data5/galaxy/project/data/total_m6a_peak"
input_dir = "/data5/galaxy/project/data/input_data/input_bed"
input_cycle_dir = "/data5/galaxy/project/data/input_data/input_cycle"
####################
reference_genome = "/data/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
motif_file = "/data/database/tf_motif/motif_matrix.txt"
####################
markov_file_dir = "/data5/galaxy/project/tf_analysis/motif_analysis/markov_bg_files"
result_dir = "/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result"
for i_dir in [markov_file_dir, result_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)
###########################################################################################


def transform_bed_to_fa(bed_file):
    fasta_file = "%s.fa" % bed_file.split(".bed")[0].split("_discrete")[0]
    os.system("bedtools getfasta -fi %s -bed %s -fo %s" % (reference_genome, bed_file, fasta_file))
    # bedtools getfasta -fi /data/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa -bed all_tissues.bed -fo all_tissues.fa


def generate_fa_file():
    ip_list = glob.glob("%s/*.bed" % ip_dir)
    print(ip_list)
    input_list = glob.glob("%s/*.bed" % input_dir)
    print(input_list)
    cycle_list = glob.glob("%s/*.bed" % input_cycle_dir)
    # print(cycle_list)
    total_bed = ip_list + input_list + cycle_list
    pool_fa = Pool()
    for bed in total_bed:
        pool_fa.apply_async(transform_bed_to_fa, (bed, ))
    pool_fa.close()
    pool_fa.join()


# Fasta-get-markov estimates a Markov model from a FASTA file of sequences. It removes ambiguous characters before
# computing the model. The model is based on both strands when using a complementable alphabet.
# Each tissue map to one markov file.

def fasta_get_markov(i_fasta, i_markov):
    os.system("fasta-get-markov %s %s" % (i_fasta, i_markov))


def get_markov_background():
    input_fa_list = glob.glob("%s/*.fa" % input_dir)
    pool_mk = Pool()
    for input_fa in input_fa_list:
        markov_file = os.path.join(markov_file_dir, "%s.txt" % os.path.basename(input_fa).split(".")[0])
        pool_mk.apply_async(fasta_get_markov, (input_fa, markov_file))
    pool_mk.close()
    pool_mk.join()


# match ip and input
# one ip map to 10 input files
def get_fa_list():
    ip_fa_list = glob.glob("%s/*.fa" % ip_dir)
    input_fa_list = glob.glob("%s/*.fa" % input_cycle_dir)
    # input_fa_list = [x for x in glob.glob("%s/*.fa" % input_cycle_dir) if int(x.split("_")[-1].split(".fa")[0]) <= 10]
    map_dict = {}
    for ip_fa in ip_fa_list:
        tissue_name = os.path.basename(ip_fa).split(".")[0].lower()
        for input_fa in input_fa_list:
            i_tissue_name = input_fa.split("-")[1].split("_")[0].lower()
            if i_tissue_name == tissue_name:
                map_dict[ip_fa] = map_dict.get(ip_fa, []) + [input_fa]
    return map_dict


def auto_decide_ip_or_input(i_fa):
    if len(os.path.basename(i_fa).split("-")) > 1:
        out_dir = "%s/%s" % (result_dir, os.path.basename(i_fa).split("-")[1].split(".")[0])
        tissue_name = os.path.basename(i_fa).split("-")[1].split("_")[0].lower()
    else:
        out_dir = "%s/ip_%s" % (result_dir, os.path.basename(i_fa).split(".")[0])
        tissue_name = os.path.basename(i_fa).split(".")[0].lower()
    return out_dir, tissue_name


def identify_motif_from_fa(i_fa, markov_files):
    out_dir, tissue_name = auto_decide_ip_or_input(i_fa)
    markov_file = [markov for markov in markov_files if tissue_name in markov][0]
    # os.system("fimo --bgfile %s --thresh 1e-5 --oc %s %s %s" % (markov_file, out_dir, motif_file, ip_fa))
    os.system("fimo --bgfile %s --oc %s %s %s" % (markov_file, out_dir, motif_file, i_fa))


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #
    generate_fa_file()
    get_markov_background()
    markov_file_list = glob.glob("%s/*.txt" % markov_file_dir)
    mapping_dict = get_fa_list()
    pool = Pool(processes=20)
    for ip_fasta, input_fasta_list in mapping_dict.items():
        # print(ip_fasta)
        # print(input_fasta_list)
        input_fasta_list.append(ip_fasta)   #
        for input_fasta in input_fasta_list:
            i_out_dir, tissue_n = auto_decide_ip_or_input(input_fasta)
            # if not os.path.exists(i_out_dir):
            pool.apply_async(identify_motif_from_fa, (input_fasta, markov_file_list))
            # identify_motif_from_fa(input_fasta, markov_file_list)
    pool.close()
    pool.join()
    #
    print(start_time)
    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(end_time)
