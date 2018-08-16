#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
os.chdir("/data4/database/tf_motif")

motif_file = "JASPAR_HOCOMOCO_HumanTF_Uniprobe_motif.txt"
result_file = "motif_matrix.txt"


result_list = []
with open(motif_file, 'r') as f:
    for line in f.readlines():
        if line.startswith("MOTIF"):
            info = line.strip().split(" ")
            new_line = "\t".join([info[0], info[-1].split("_")[0].upper()])
            result_list.append(new_line)
        else:
            result_list.append(line.strip())

# read into final_dict
final_list, final_dict = [], {}
i = 0
# read header
while not result_list[i].startswith("MOTIF"):
    final_list.append(result_list[i])
    i += 1
# read motif and matrix
while i < len(result_list):
    if result_list[i].startswith("MOTIF"):
        # only keep TF name, so if have repeat key, the later will replace the header.
        name = result_list[i].split("\t")[-1].upper()
        final_dict[name] = []
        i += 1
        while i < len(result_list) and not result_list[i].startswith("MOTIF"):
            final_dict[name].append(result_list[i].strip())
            # print(result_list[i])
            i += 1
    else:
        i += 1

for title, values in final_dict.items():
    print(title)
    final_list.append("%s\t%s" % ("MOTIF", title.strip()))
    for x in values:
        final_list.append(x.strip())
with open(result_file, 'w') as fw:
    for line in final_list:
        fw.write(line.strip() + "\n")
