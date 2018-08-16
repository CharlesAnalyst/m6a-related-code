#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd

#############################################################################
file_disease_dir = "/data5/galaxy/project/tf_pathway/TF_disease/file_result"
file_list = glob.glob("%s/*.txt" % file_disease_dir)
snp_disease_dir = "/data5/galaxy/project/tf_pathway/TF_disease/snp_result"
snp_list = glob.glob("%s/*.txt" % snp_disease_dir)
result_dir = "/data5/galaxy/project/tf_pathway/TF_disease/compare_disease"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#####################################################


def get_snp_by_file(i_file):
    snp_file = ""
    tissue = os.path.basename(i_file).split(".txt")[0]
    for snp in snp_list:
        if tissue in snp:
            snp_file = snp
            break
    return snp_file


def generate_dict_from_file(in_file):
    result_dict = {}
    df = pd.read_table(in_file, sep="\t", header=None, index_col=0, names=["disease"])
    for name, values in df.iterrows():
        result_dict[name] = result_dict.get(name, []) + [values["disease"]]
    return result_dict


def main_method():
    title = "transcript factor\tdisease(file)\tdisease(snp)\tcommon disease\tcommon disease number\n"
    for i_f in file_list:
        result_file = os.path.join(result_dir, os.path.basename(i_f))
        with open(result_file, 'w') as fw:
            fw.write(title)
            i_snp = get_snp_by_file(i_f)
            file_dict = generate_dict_from_file(i_f)
            snp_dict = generate_dict_from_file(i_snp)
            key_list = list(set(list(file_dict.keys()) + list(snp_dict.keys())))
            for i_key in key_list:
                if i_key in file_dict and i_key in snp_dict:
                    values = file_dict[i_key]
                    snp_dis_list = snp_dict[i_key]
                    set_a, set_b = set(values), set(snp_dis_list)
                    overlap_list = list(set_a & set_b)
                    result_line = "%s\t%s\t%s\t%s\t%d\n" % (i_key, ";".join(values), ";".join(snp_dis_list), ";".join(overlap_list), len(overlap_list))
                elif i_key in file_dict:
                    values = file_dict[i_key]
                    result_line = "%s\t%s\t%s\t%s\t%d\n" % (i_key, ";".join(values), "", "", 0)
                else:
                    snp_dis_list = snp_dict[i_key]
                    result_line = "%s\t%s\t%s\t%s\t%d\n" % (i_key, "", ";".join(snp_dis_list), "", 0)
                fw.write(result_line)


if __name__ == '__main__':
    main_method()
