#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import pandas as pd

##############################################################################################################
id_transform_file = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/04_intersect/id_trans.txt"
##############################################################


def main_method():
    work_dir = "/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result"
    result_file = "%s/zscore_filter-by-RPKM_result.txt" % work_dir

    info_table = "/data5/xlj/new_tissue_m6A/Macs2/Filter/Merge-Dup/total_info.txt"  # "entrez_id"
    df = pd.read_table(info_table, sep="\t")
    col_list = ["entrez_id"] + [col for col in df.columns if "_exp" in col]
    df_sub = df[col_list]   # expression
    # get trans_id by total express file(contain total gene)
    trans_id = pd.read_table(id_transform_file, sep="\t")

    motif_file = "%s/zscore_filter-by-cutoff.txt" % work_dir
    df_motif = pd.read_table(motif_file, sep="\t")

    df_merge = pd.merge(df_motif, trans_id, left_on="motif_id", right_on="SYMBOL", how="left")
    df_merge = df_merge.dropna(how="any")
    print(df_merge)
    del df_merge["motif_id"]
    del df_merge["SYMBOL"]
    #
    df_expre = df_sub.set_index("entrez_id")
    df_motif = df_merge.set_index("ENSEMBL")
    #
    df_expre = df_expre.sort_index(axis=1)
    df_motif = df_motif.sort_index(axis=1)
    # print(df_motif)
    #
    col_list = df_motif.columns
    index_list = df_expre.index
    with open(result_file, 'w') as fw:
        fw.write("ensembl_id\t%s\n" % "\t".join(col_list))
        for name, values in df_motif.iterrows():
            # print("name: %s" % name)
            if name in index_list:
                values_expre = df_expre.loc[name, :]
                line_list = []
                for i in range(len(col_list)):
                    # if values[i] > 1.1065 and values_expre[i] > values_expre.mean():
                    if values_expre[i] > values_expre.mean():
                        line_list.append(values[i])
                    else:
                        line_list.append(0.0)
                fw.write("%s\t%s\n" % (name, "\t".join([str(x) for x in line_list])))
            else:
                print("non: %s" % name)


main_method()
