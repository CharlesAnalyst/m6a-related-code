#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd

###################################################################################
tf_dir = "/data5/galaxy/project/tf_pathway/m6aSNP_TF_intersect/TF"
tf_file_list = glob.glob("%s/*.txt" % tf_dir)
tf2dis_file = "/data4/database/DisGeNET_db/curated_gene_disease_associations.tsv"
#########################
tf_snp_dir = "/data5/galaxy/project/tf_pathway/m6aSNP_TF_intersect/tf_snp"
tf_snp_file_list = glob.glob("%s/*.txt" % tf_snp_dir)
snp2dis_file = "/data4/database/DisGeNET_db/curated_variant_disease_associations.tsv"
########################
result_dir = "/data5/galaxy/project/tf_pathway/TF_disease"
file_dir, snp_dir = os.path.join(result_dir, "file_result"), os.path.join(result_dir, "snp_result")
for i_dir in [result_dir, file_dir, snp_dir]:
    if not os.path.exists(i_dir):
        os.makedirs(i_dir)
###################################################################################


def get_disease_from_file():
    tf_dis_dict = {}
    df_gene2disease = pd.read_table(tf2dis_file, sep="\t")
    for name, values in df_gene2disease.iterrows():
        i_key, i_value = values.geneSymbol, values.diseaseName
        tf_dis_dict[i_key] = tf_dis_dict.get(i_key, []) + [i_value]
    for tf_file in tf_file_list:
        out_file = os.path.join(file_dir, os.path.basename(tf_file))
        if os.path.exists(out_file):
            os.remove(out_file)
        with open(out_file, 'a+') as fw:
            tf_list = pd.read_table(tf_file, sep="\t").iloc[:, 0].tolist()
            for tf in tf_list:
                if tf.upper() in tf_dis_dict:
                    disease_list = tf_dis_dict[tf.upper()]
                    for disease in disease_list:
                        fw.write("%s\t%s\n" % (tf, disease))
                else:
                    print("%s not in File!" % tf)


def get_disease_from_overlapped_snp():
    snp_dis_dict = {}
    df_snp2disease = pd.read_table(snp2dis_file, sep="\t")
    for name, values in df_snp2disease.iterrows():
        i_key, i_value = values.snpId, values.diseaseName
        snp_dis_dict[i_key] = snp_dis_dict.get(i_key, []) + [i_value]
    for tf_snp in tf_snp_file_list:
        out_file = os.path.join(snp_dir, os.path.basename(tf_snp))
        if os.path.exists(out_file):
            os.remove(out_file)
        tf_dis_dict = {}
        df = pd.read_table(tf_snp, sep="\t", header=None, names=["tf", "snp"], index_col=0)
        for name, values in df.iterrows():
            if values["snp"] in snp_dis_dict:
                disease = snp_dis_dict[values["snp"]]
                tf_dis_dict[name] = tf_dis_dict.get(name, []) + disease
        with open(out_file, 'w') as fw:
            for tf, disease_list in tf_dis_dict.items():
                uniq_list = list(set(disease_list))
                for disease in uniq_list:
                    fw.write("%s\t%s\n" % (tf, disease))


get_disease_from_file()
get_disease_from_overlapped_snp()
