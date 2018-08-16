#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import subprocess

########################################################################################################################
gwas_snp = "/data5/galaxy/project/GWAS_db/snp_classification/merge/total_gwas_GRCh38.bed"
eQTL_dir = "/data5/galaxy/project/snp_analysis/snp_db/GTEx_Analysis_v7_eQTL/format_bed/GRCh37_to_GRCh38"
eQTL_list = glob.glob("%s/*.bed" % eQTL_dir)
m6a_dir = "/data5/galaxy/project/data/total_m6a_peak"
m6a_list = glob.glob("%s/*.bed" % m6a_dir)
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/create_big-table"
tissue_list = ["BRAIN", "HEART", "LIVER", "LUNG", "MUSCLE", "STOMACH"]
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#################################################################################################################


def get_each_tissue_files(file_list):
    tissue_file_dict = {}
    for i_file in file_list:
        tissue = os.path.basename(i_file).split(".bed")[0].split("_")[0].upper()
        tissue_file_dict[tissue] = tissue_file_dict.get(tissue, []) + [i_file]
    return tissue_file_dict


def intersect_with_m6a(in_bed, m6a_bed):
    result_dict = {}
    awk_command = 'print $4 "\t" $5 ":" $6 "-" $7'
    command = "bedtools intersect -a %s -b %s -wa -wb| awk '{%s}'" % (in_bed, m6a_bed, awk_command)
    sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    results = sub_p.communicate()[0].split("\n")
    for line in results:
        if not line == "":
            # print(line)
            result_dict[str(line).split("\t")[0]] = str(line).split("\t")[1]
    return result_dict


def main_method():
    snp_dict, m6a_dict = get_each_tissue_files(eQTL_list), get_each_tissue_files(m6a_list)
    for tissue in tissue_list:
        print(tissue)
        snp_list, m6a_bed = snp_dict[tissue], m6a_dict[tissue][0]
        for snp_bed in snp_list:
            print(snp_bed, m6a_bed)
            result_file = os.path.join(result_dir, os.path.basename(snp_bed).replace(".bed", ".txt"))
            eQTL_dict, gwas_dict = intersect_with_m6a(snp_bed, m6a_bed), intersect_with_m6a(gwas_snp, m6a_bed)
            total_snp = list(set(list(eQTL_dict.keys()) + list(gwas_dict.keys())))
            with open(result_file, 'w') as fw:
                fw.write("rs_id\tLocation\teQTL_SNP\tGWAS_SNP\tBoth\n")
                for snp in total_snp:
                    result_line = [snp, "False", "False", "False", "False"]
                    if snp in eQTL_dict:
                        result_line[1] = eQTL_dict[snp]
                        result_line[2] = "True"
                        if snp in gwas_dict:
                            result_line[-2] = "True"
                            result_line[-1] = "True"
                    elif snp in gwas_dict:
                        result_line[1] = gwas_dict[snp]
                        result_line[-2] = "True"
                    fw.write("\t".join(result_line) + "\n")


if __name__ == '__main__':
    main_method()
