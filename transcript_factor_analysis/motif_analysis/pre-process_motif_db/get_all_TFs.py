#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


in_file = "D:\\Project\\Motif\\new_download_data\\MEME_HOCOMOCO\\HOCOMOCOv10_HUMAN_mono_meme_format.meme"
out_file = "D:\\Project\\Motif\\new_download_data\\MEME_HOCOMOCO\\HOCOMOCOv10_motif_names.txt"
a = 'D:\software\复制文件路径B2.exe'

with open(out_file, 'w') as fw:
    fw.write("Motif\tmotif_name\n")
    with open(in_file, 'r', encoding="utf-8") as f:
        for line in f.readlines():
            if line.startswith("MOTIF"):
                fw.write(line.strip().replace(" ", "\t") + "\n")


