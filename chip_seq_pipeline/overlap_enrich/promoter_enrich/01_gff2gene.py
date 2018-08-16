#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import pandas as pd

##############################
gff, genome_region_bed = "", ""
species = "human"
if species == "mouse":
    gff = "/data/database/GRCm38/GENCODE/gencode.vM15.annotation.gff3"
    genome_region_bed = "/data/database/GRCm38/GENCODE/Genes_ensembl.bed"
elif species == "human":
    gff = "/data/database/GRCh38/GENCODE/gencode.v27.annotation.gff3"
    genome_region_bed = "/data/database/GRCh38/GENCODE/Genes_ensembl_dot.bed"
###############################


def get_genome_region():
    title = ["chr", "HAVANA", "type", "start", "end", "a", "strand", "b", "others"]
    data = pd.read_table(gff, sep="\t", skiprows=7, names=title)
    sub_data = data[data["type"] == "gene"]
    id_list = [i.split("gene_id=")[1].split(";")[0] for i in list(sub_data["others"])]
    sub_data.index = id_list
    sub_data.reset_index(inplace=True)
    result_data = sub_data[["chr", "start", "end", "index", "a", "strand"]]
    result_data[["start", "end"]] = result_data[["start", "end"]].astype("int")
    result_data.to_csv(genome_region_bed, sep="\t", index=False, header=False)


if __name__ == '__main__':
    get_genome_region()

"""
def reformat():
    with open(transcript_bed, 'r') as f:
        with open(reformat_file, 'w') as fw:
            for line in f.readlines():
                info = line.split("\t")
                chr_id, start, end, name = info[0], info[1], info[2], info[3]
                fw.write("\t".join([name.split(".")[0], chr_id, start, end]) + "\n")

reformat()
"""