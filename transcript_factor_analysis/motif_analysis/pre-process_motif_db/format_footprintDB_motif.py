#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob

jaspar_dir = "/data4/galaxy/project/motif_match/database/new_download_data/JASPAR/DBM"
uniprobe_dir = "/data4/galaxy/project/motif_match/database/new_download_data/uniprob/DBM"
humanTF_dir = "/data4/galaxy/project/motif_match/database/new_download_data/HumanTF1/DBM"
total_jaspar_file = "/data4/galaxy/project/motif_match/database/new_download_data/JASPAR/total_motifs.txt"
total_uniprobe_file = "/data4/galaxy/project/motif_match/database/new_download_data/uniprob/total_motifs.txt"
total_humanTF_file = "/data4/galaxy/project/motif_match/database/new_download_data/HumanTF1/total_motifs.txt"


def generate_file_header():
    result_list = ["MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 "
                   "G 0.25 T 0.25\n\n"]
    return result_list


def format_format_to_meme(in_file):
    result_list, first_line, second_line, matrix_list = [], "", "", []
    with open(in_file, 'r') as f:
        motif_name = f.readline().split("|")[1].split("(")[0]
        contents = f.readlines()
        for i in range(len(contents)):
            if contents[i].startswith("NA"):
                alternate_name = contents[i].split("NA")[1].strip()
                first_line = "MOTIF %s %s" % (motif_name, alternate_name)
            if contents[i].startswith("P0"):
                i += 1
                while not contents[i].startswith("XX"):
                    val_list = [float(val) for val in contents[i].split(" ")[1:-1] if val != ""]
                    new_val_list = [str(val/sum(val_list)) for val in val_list]
                    matrix_list.append(" ".join(new_val_list))
                    i += 1
                second_line = "letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0.0" % len(matrix_list)
        result_list.append(first_line)
        result_list.append(second_line)
        for line in matrix_list:
            result_list.append(line)
        result_list.append("\n")
    return result_list


def write_to_file(out_file, in_list):
    with open(out_file, 'w') as fw:
        for line in in_list:
            fw.write(line + "\n")


def transform_and_merge_motifs(file_dir, out_file):
    os.chdir(file_dir)
    result_list = []
    motif_files = glob.glob("*.motif.transfac")
    header_list = generate_file_header()
    result_list += header_list
    for file in motif_files:
        matrix_list = format_format_to_meme(file)
        result_list += matrix_list
    write_to_file(out_file, result_list)


if __name__ == "__main__":
    # transform_and_merge_motifs(jaspar_dir, total_jaspar_file)
    # transform_and_merge_motifs(uniprobe_dir, total_uniprobe_file)
    transform_and_merge_motifs(humanTF_dir, total_humanTF_file)

