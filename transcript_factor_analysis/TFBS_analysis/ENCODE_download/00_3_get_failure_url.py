#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os

os.chdir("/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/raw_data")
os.system("md5sum --check md5_sum.txt > check.log")
log_file = "/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/raw_data/check.log"
power_shell = "/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/power_shell.txt"
result_file = "/data5/galaxy/project/TF_narrowPeak/raw_data/total_TF_narrowPeak/power_shell_second.txt"

failure_file_list = []
with open(log_file, 'r') as f:
    for line in f.readlines():
        if "OK" not in line:
            failure_file_list.append(line.split(":")[0])

with open(result_file, 'w') as fw:
    with open(power_shell, 'r') as f:
        title = f.readline()
        contents = f.readlines()
        fw.write(title)
        for line in contents:
            if line.split("/")[-1].split("',")[0] in failure_file_list:
                fw.write(line)
