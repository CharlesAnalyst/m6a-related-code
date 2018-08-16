#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
from collections import Counter

############################################################
fimo_dir = "/data5/galaxy/project/tf_pathway/fimo_result"
result_dir = "/data5/galaxy/project/tf_pathway/m6a_snp_TF"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
###########################################################

fimo_list = glob.glob("%s/*/fimo.txt" % fimo_dir)
for fimo in fimo_list:
    out_file = fimo.replace("fimo.txt", "TF_count.txt")
    df = pd.read_table(fimo)
    tf_list = df["# motif_id"].tolist()
    count = Counter(tf_list)
    sorted_count = sorted(count.items(), key=lambda a: a[1], reverse=True)
    with open(out_file, 'w') as fw:
        for x, c in sorted_count:
            fw.write(x + "\t" + str(c) + "\n")
