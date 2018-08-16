#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import re
import math
import subprocess
import pandas as pd


##########################################################################################
CUTOFF = 0.35
UPSTREAM = 2000
DOWNSTREAM = 2000
PROMOTER_UPSTREAM = 2000
#######################
reference_genome = "/data/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
#
promoter_bed = "/data5/galaxy/project/data/promoter/human/2k_100/genes_promoters.bed"
mettl3_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"
#
result_dir = "/data5/galaxy/project/mettl3_enrich/mettl3_enrich_high_CpG"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#######################
# ip_bw = "/data5/galaxy/project/mettl3_enrich/clean_bam/plot_profile/ip_coverage.bw"
# input_bw = "/data5/galaxy/project/mettl3_enrich/clean_bam/plot_profile/input_coverage.bw"
##########################################################################################


def calculate_CpG_density_for_all():
    print("enrich CpG in promoter")
    title_list, score_list = enrich_CpG_in_promoter()
    print("cluster promoter into")
    df_pro_high, df_pro_low = cluster_promoter_into_high_low(title_list, score_list)
    print("get tss region")
    region_bed_high, region_bed_low = get_tss_region(df_pro_high, "high"), get_tss_region(df_pro_low, "low")
    print("plot")
    # plot_profile(region_bed_high, region_bed_low)   # , plot_profile(region_bed_low, "low")


def enrich_CpG_in_promoter():
    str_seq_list, title_list, score_list = get_sequence_from_bed().split("\n"), [], []
    for i in range(0, len(str_seq_list)-1, 2):
        title, seq = str_seq_list[i], str_seq_list[i+1]
        title_list.append(title)
        score_list.append(calculate_CpG_density(seq))
    return title_list, score_list


def get_sequence_from_bed():
    command = "bedtools getfasta -fi %s -bed %s" % (reference_genome, promoter_bed)
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    str_sequence = sub_p.communicate()[0]
    return str_sequence


def calculate_CpG_density(sequence):
    seq = sequence.lower()
    cg_num, c_num, g_num = len(re.findall("cg", seq)), len(re.findall("c", seq)), len(re.findall("g", seq))
    score = cg_num / (math.pow((c_num + g_num) / 2.0, 2) / (len(seq) * 1.0))
    return score


def cluster_promoter_into_high_low(title_list, score_list):
    # high_index, low_index = calculate_kmeans.k_means(2, score_list)
    high_index = [i for i in range(len(score_list)) if score_list[i] > CUTOFF]
    low_index = [i for i in range(len(score_list)) if score_list[i] <= CUTOFF]
    print(len(high_index), len(low_index))
    high_title, low_title = [title_list[x] for x in high_index], [title_list[x] for x in low_index]
    df_pro_high, df_pro_low = get_promoter_according_title(high_title), get_promoter_according_title(low_title)
    return df_pro_high, df_pro_low


def get_promoter_according_title(title_list):
    chro_list = [title.split(":")[0].split(">")[1] for title in title_list]
    start_list = [int(title.split(":")[1].split("-")[0]) for title in title_list]
    end_list = [int(title.strip().split("-")[-1]) for title in title_list]
    df_title = pd.DataFrame({"chr": chro_list, "start": start_list, "end": end_list}, columns=["chr", "start", "end"])
    df = pd.read_table(promoter_bed, sep="\t", header=None, names=["chr", "start", "end", "d", "e", "strand"])
    df_select = pd.merge(df_title, df, how="left", on=["chr", "start", "end"]).dropna(how="any")
    df_select = df_select[["chr", "start", "end", "d", "e", "strand"]]
    return df_select


def get_tss_region(df, data_type):
    df_pos, df_neg = df[df.loc[:, "strand"] == "+"], df[df.loc[:, "strand"] == "-"]
    #
    df_neg["tss_s"], df_neg["tss_e"] = ((df_neg["end"] - PROMOTER_UPSTREAM) - DOWNSTREAM), ((df_neg["end"] - PROMOTER_UPSTREAM) + UPSTREAM)
    df_neg_bed = df_neg[["chr", "tss_s", "tss_e", "d", "e", "strand"]]
    #
    df_pos["tss_s"], df_pos["tss_e"] = ((df_pos["start"] + PROMOTER_UPSTREAM) - UPSTREAM), ((df_pos["start"] + PROMOTER_UPSTREAM) + DOWNSTREAM)
    df_pos_bed = df_pos[["chr", "tss_s", "tss_e", "d", "e", "strand"]]
    #
    df_tss_region = pd.concat([df_pos_bed, df_neg_bed]).sort_values(["chr", "tss_s"])
    df_tss_region = df_tss_region[df_tss_region["tss_s"] > 0]
    tss_region_bed = os.path.join(result_dir, "tss_%s.bed" % data_type)
    df_tss_region.to_csv(tss_region_bed, sep="\t", header=None, index=False)
    # tss_string = str(df_tss_region.to_string(header=False, index=False, col_space=4))
    # tss_region_string = "\n".join(["\t".join(x.split()) for x in tss_string.split("\n")])
    return tss_region_bed


"""
def plot_profile(high_region_bed, low_region_bed):
    result_file = os.path.join(result_dir, "mettl3_tss.gz")
    command = "computeMatrix reference-point --referencePoint center --numberOfProcessors max/2 -S %s -R %s %s -a %s -b %s -o %s" % (ip_bw, high_region_bed, low_region_bed, UPSTREAM, DOWNSTREAM, result_file)
    os.system(command)
    figure_file = os.path.join(result_dir, os.path.basename(result_file).replace(".gz", ".pdf"))
    command = "plotProfile -m %s -out %s --plotTitle ''" % (result_file, figure_file)
    os.system(command)
"""


if __name__ == '__main__':
    calculate_CpG_density_for_all()
