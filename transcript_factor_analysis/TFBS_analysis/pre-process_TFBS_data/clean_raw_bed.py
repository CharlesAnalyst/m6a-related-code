#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import pandas as pd

###################################################################
file_dir = "/data3/xs/tissue_m6a/2018.1/TF/ENCODE_PEAKS/"
bed_list = glob.glob("%s/*/*.bed" % file_dir)
result_dir = "/data5/galaxy/project/tf_analysis/TFBS_analysis/data"
###################################################################


def transcript_factor_from_bed(in_bed):
    tf_name = in_bed.split("/")[-2].upper()
    return tf_name


def get_trans_factor_dict():
    transcript_factor_dict = {}
    for bed in bed_list:
        tf_name = transcript_factor_from_bed(bed)
        transcript_factor_dict[tf_name] = transcript_factor_dict.get(tf_name, []) + [bed]
    return transcript_factor_dict


def process_one_directory(i_bed_list, out_dir):
    chr_list = ["chr%d" % i for i in range(1, 23)]
    chr_list.append("chrX")
    chr_list.append("chrY")
    # extract only three columns
    for file in i_bed_list:
        print(file)
        df = pd.read_table(file, sep="\t", header=None)
        df_sub = df.iloc[:, :3]
        #
        result_list = []
        for name, values in df_sub.iterrows():
            # extract chr in chr_list[1-22,X,Y]
            if values[0] in chr_list:
                result_list.append(values.tolist())
        df_result = pd.DataFrame(result_list, columns=["chr", "start", "end"])
        out_file = os.path.join(out_dir, os.path.basename(file))
        df_result.to_csv(out_file, sep="\t", header=None, index=None)


if __name__ == "__main__":
    tf_dict = get_trans_factor_dict()
    for tf, beds in tf_dict.items():
        i_out_dir = os.path.join(result_dir, tf)
        if not os.path.exists(i_out_dir):
            os.makedirs(i_out_dir)
        process_one_directory(beds, i_out_dir)
