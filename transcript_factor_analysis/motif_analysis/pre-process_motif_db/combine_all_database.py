#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os

human_tf = "/data4/galaxy/project/motif_match/database/new_download_data/HumanTF1/total_motifs.txt"
jaspar = "/data4/galaxy/project/motif_match/database/new_download_data/JASPAR/total_motifs.txt"
hocomoco = "/data4/galaxy/project/motif_match/database/new_download_data/MEME_HOCOMOCO" \
           "/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
uniprobe = "/data4/galaxy/project/motif_match/database/new_download_data/uniprob/total_motifs.txt"


for file in [human_tf, jaspar, hocomoco, uniprobe]:
    out_file = os.path.join(os.path.split(file)[0], "motif_names.txt")
    with open(out_file, 'w') as fw:
        with open(file, 'r') as f:
            for line in f.readlines():
                if line.startswith("MOTIF"):
                    fw.write("%s\n" % line.split(" ")[1].strip())

