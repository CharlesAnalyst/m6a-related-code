#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
import subprocess


######################################################################################################
data_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/gain_or_loss/01_predict_result_hg19"
info_list = glob.glob("%s/*.txt" % data_dir)
postfix = ".txt"
m6a_dir = "/data5/galaxy/project/data/total_m6a_peak/hg19_version"
m6a_list = glob.glob("%s/*.bed" % m6a_dir)
tissue_list = ["BRAIN", "HEART", "LIVER", "LUNG", "MUSCLE", "STOMACH"]
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/gain_or_loss/02_parse_result"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#######################################################################################################
raw_gain_list, raw_loss_list, alt_gain_list, alt_loss_list = [], [], [], []
raw_gain_file, raw_loss_file, alt_gain_file, alt_loss_file = os.path.join(result_dir, "raw_gain.fa"), \
                                                             os.path.join(result_dir, "raw_loss.fa"), \
                                                             os.path.join(result_dir, "alt_gain.fa"), \
                                                             os.path.join(result_dir, "alt_loss.fa")
result_file_list = [raw_gain_file, raw_loss_file, alt_gain_file, alt_loss_file]


def main_method():
    info_dict, m6a_dict = get_each_tissue_files(info_list), get_each_tissue_files(m6a_list)
    for tissue in tissue_list:
        print(tissue)
        info_file_list, m6a_bed = info_dict[tissue], m6a_dict[tissue][0]
        for info_bed in info_file_list:
            print(info_bed, m6a_bed)
            prefix = os.path.basename(info_bed).split(".txt")[0]
            df = pd.read_table(info_bed, sep="\t", index_col=False)
            df_gain = df[df["Mutation event"] == "Functional Gain"]
            df_loss = df[df["Mutation event"] == "Functional Loss"]
            global raw_gain_list, raw_loss_list, alt_gain_list, alt_loss_list
            ref_seq_list, alt_seq_list = process_df(df_gain, m6a_bed, prefix)
            raw_gain_list += ref_seq_list
            alt_gain_list += alt_seq_list
            ref_seq_list, alt_seq_list = process_df(df_loss, m6a_bed, prefix)
            raw_loss_list += ref_seq_list
            alt_loss_list += alt_seq_list
    result_list = [raw_gain_list, raw_loss_list, alt_gain_list, alt_loss_list]
    write_to_file(result_list)


def get_each_tissue_files(file_list):
    tissue_file_dict = {}
    for i_file in file_list:
        tissue = os.path.basename(i_file).split(".bed")[0].split(".txt")[0].split("_")[0].upper()
        tissue_file_dict[tissue] = tissue_file_dict.get(tissue, []) + [i_file]
    return tissue_file_dict


def process_df(df_in, m6a_bed, prefix):
    ref_seq_list, alt_seq_list = [], []
    # df_in = df_in[["Gene symbol", "Chromosome", "Position", "Strand", "Related mutations"]]
    df_in["rs_ID"] = df_in["Related mutations"].str.split("|").str[2]
    df_in["name"] = df_in["rs_ID"] + ";" + df_in["Reference sequence"] + ";" + df_in["Mutation sequence"]
    df_in["end"] = df_in["Position"].astype(int) + 1
    df_bed = df_in[["Chromosome", "Position", "end", "name"]].sort_values(["Chromosome", "Position"])
    df_list = []
    for name, values in df_bed.iterrows():
        df_list.append("\t".join([str(x) for x in values]))
    content = ("\n".join(df_list)).encode("utf-8")
    # print(content)
    command = "bedtools intersect -a stdin -b %s -wa | sort -k1,1 -k2,2n | uniq" % m6a_bed
    sub_p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    m6a_seq_list = str(sub_p.communicate(content)[0].decode("utf-8")).split("\n")
    print(m6a_seq_list)
    m6a_seq_list.remove("")
    if len(m6a_seq_list) > 0:
        ref_seq_list = [">%s_%d\n%s\n" % (prefix, i, m6a_seq_list[i].split()[3].split(";")[1]) for i in range(len(m6a_seq_list))]
        alt_seq_list = [">%s_%d\n%s\n" % (prefix, i, m6a_seq_list[i].split()[3].split(";")[-1]) for i in range(len(m6a_seq_list))]
    return ref_seq_list, alt_seq_list


def write_to_file(result_list):
    for i in range(len(result_file_list)):
        with open(result_file_list[i], 'w') as fw:
            fw.writelines(result_list[i])


if __name__ == '__main__':
    main_method()
