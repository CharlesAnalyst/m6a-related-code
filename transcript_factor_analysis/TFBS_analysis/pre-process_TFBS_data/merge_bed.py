#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob


########################################################################
chip_peak_dir = "/data5/galaxy/project/tf_analysis/TFBS_analysis/data"
bed_list = glob.glob("%s/*/*.bed" % chip_peak_dir)
result_dir = "/data5/galaxy/project/tf_analysis/TFBS_analysis/ip_merge"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#########################################################################


def transcript_factor_from_bed(in_bed):
    tf_name = in_bed.split("/")[-2].upper()
    return tf_name


def get_trans_factor_dict():
    transcript_factor_dict = {}
    for bed in bed_list:
        tf_name = transcript_factor_from_bed(bed)
        transcript_factor_dict[tf_name] = transcript_factor_dict.get(tf_name, []) + [bed]
    return transcript_factor_dict


def get_merged_bed(transcript_factor_dict):
    for transcript_factor, beds in transcript_factor_dict.items():
        print(transcript_factor)
        merged_bed = os.path.join(result_dir, "%s.bed" % transcript_factor)
        os.system("cat %s | sort -k1,1 -k2,2n - | bedtools merge -i - > %s" % (" ".join(beds), merged_bed))
        transcript_factor_dict[transcript_factor] = merged_bed


if __name__ == "__main__":
    tf_dict = get_trans_factor_dict()
    get_merged_bed(tf_dict)
