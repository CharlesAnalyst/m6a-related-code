#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob

bam_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/unique_bam/"
postfix = "_sorted_unique.bam"
result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/macs2_peak"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def main():
    result_dict = generate_ip_control_pair()
    for name, file_dict in result_dict.items():
        ip_bam, input_bam = file_dict["ip"], file_dict["input"]
        print(ip_bam)
        print(input_bam)
        call_peak(ip_bam, input_bam, name)


def generate_ip_control_pair():
    type_dict, result_dict = {}, {}
    bam_list = glob.glob("%s/*%s" % (bam_dir, postfix))
    for bam in bam_list:
        prefix = os.path.basename(bam).split("_")[0]
        type_dict[prefix] = type_dict.get(prefix, []) + [bam]
    for prefix, file_list in type_dict.items():
        result_dict[prefix] = {}
        for i_file in file_list:
            if "_IP_" in i_file:
                result_dict[prefix]["ip"] = i_file
            elif "_In_" not in i_file:
                result_dict[prefix]["input"] = i_file
    return result_dict


def call_peak(ip_bam, input_bam, exp_name):
    command = "macs2 callpeak -B -t %s -c %s --outdir %s -q 0.05 --nomodel -n %s -g mm" % (ip_bam, input_bam, result_dir, exp_name)
    os.system(command)
    in_bed, result_bed = os.path.join(result_dir, "%s_peaks.narrowPeak" % exp_name), os.path.join(result_dir, "%s.bed" % exp_name)
    os.system("cut -f 1-3 %s > %s" % (in_bed, result_bed))


if __name__ == '__main__':
    main()