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
species = "hg38"
#######################
reference_genome = "/data/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
promoter_bed = "/data5/galaxy/project/data/promoter/human/3k_3k/gene_promoter.bed"
gene_bed = "/data/database/GRCh38/GENCODE/Genes_ensembl.bed"
result_dir = "/data5/galaxy/project/CpG_m6a_motif/fasta_seq"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
##########################################################################################


def calculate_CpG_density_for_all():
    print("enrich CpG in promoter")
    title_list, score_list = enrich_CpG_in_promoter()
    print("cluster promoter into")
    class_promoter_into_high_low(title_list, score_list)
    # high_genes, low_genes = class_promoter_into_high_low(title_list, score_list)
    # get_gene_sequence(high_genes, os.path.join(result_dir, "high_genes.bed"))
    # get_gene_sequence(low_genes, os.path.join(result_dir, "low_genes.bed"))


def enrich_CpG_in_promoter():
    str_seq_list, title_list, score_list = str(get_sequence_from_bed()).split("\\n"), [], []
    print(len(str_seq_list))
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


def class_promoter_into_high_low(title_list, score_list):
    # high_index, low_index = calculate_kmeans.k_means(2, score_list)
    high_index = [i for i in range(len(score_list)) if score_list[i] > CUTOFF]
    low_index = [i for i in range(len(score_list)) if score_list[i] <= CUTOFF]
    print(len(high_index), len(low_index))
    high_title, low_title = ["%s;%f" % (title_list[x], score_list[x]) for x in high_index], ["%s;%f" % (title_list[x], score_list[x]) for x in low_index]
    high_pro_bed, low_pro_bed = os.path.join(result_dir, "high_CpG.bed"), os.path.join(result_dir, "low_CpG.bed")
    get_promoter_according_title(high_title, high_pro_bed)
    get_promoter_according_title(low_title, low_pro_bed)
    # high_string, low_string = homer_annotate_peak(high_pro_bed), homer_annotate_peak(low_pro_bed)
    # high_genes, low_genes = parse_annotated_region(high_string), parse_annotated_region(low_string)
    # return high_genes, low_genes


def get_promoter_according_title(title_list, result_bed):
    chro_list = [title.split(":")[0].split(">")[1] for title in title_list]
    start_list = [int(title.split(":")[1].split("-")[0]) for title in title_list]
    end_list = [int(title.strip().split("-")[-1].split(";")[0]) for title in title_list]
    score_list = [float(title.strip().split(";")[-1]) for title in title_list]
    df_title = pd.DataFrame({"chr": chro_list, "start": start_list, "end": end_list, "score": score_list})
    df = pd.read_table(promoter_bed, sep="\t", header=None, names=["chr", "start", "end", "d", "e", "strand"])
    df_select = pd.merge(df_title, df, how="left", on=["chr", "start", "end"]).dropna(how="any")
    df_select = df_select[["chr", "start", "end", "d", "score", "strand"]]
    df_select.to_csv(result_bed, sep="\t", header=False, index=False)


def homer_annotate_peak(peak_bed):
    command = "annotatePeaks.pl %s %s" % (peak_bed, species)
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    result_string = sub_p.communicate()[0]
    return result_string


def parse_annotated_region(result_string):
    line_list = str(result_string).strip().split("\n")
    col_names, content = line_list[0].split("\t"), line_list[1:]
    content_list = [line.split("\t") for line in content]
    df_content = pd.DataFrame(content_list, columns=col_names)
    gene_list = df_content["Nearest Ensembl"].tolist()
    print(len(gene_list))
    return gene_list


def get_gene_sequence(gene_list, result_bed):
    df = pd.read_table(gene_bed, sep="\t", header=None, names=["a", "b", "c", "name", "e", "f"])
    df_bed = df[df["name"].isin(gene_list)]
    df_bed.to_csv(result_bed, sep="\t", header=False, index=False)
    result_fa = result_bed.replace(".bed", ".fa")
    os.system("bedtools getfasta -fi %s -bed %s -fo %s" % (reference_genome, result_bed, result_fa))


if __name__ == '__main__':
    calculate_CpG_density_for_all()
