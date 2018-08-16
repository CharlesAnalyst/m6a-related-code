#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


##################################################################################
def figure_2_a_b():
    os.chdir("/data5/galaxy/project/trend_test")
    file_list = glob.glob("*.txt")
    for i_file in file_list:
        print(i_file)
        df = pd.read_table(i_file, sep="\t", header=None, names=["Relative expression dynamics of gene", "Raw_category"])
        df["Category"] = df["Raw_category"].str.split("-").str[1]
        sns.set_style("white")
        f, axes = plt.subplots(figsize=(10, 6))
        # sns.palplot(sns.hls_palette(8, l=1, s=.7))
        # current_palette = sns.color_palette("pastel")
        # sns.palplot(current_palette)
        sns_plot = sns.violinplot(x=df["Category"], y=df["Relative expression dynamics of gene"], palette=sns.husl_palette(8, l=.8, s=1), width=.9).get_figure()    #  scale="area", width=1.0
        axes.tick_params(direction='out', labelsize=20)
        axes.set_xlabel("", weight='bold')
        axes.set_ylabel("")
        # axes.set_ylabel("Relative expression\ndynamics of gene", fontsize=24, weight='bold', labelpad=12)
        # 去掉多余的线
        sns.despine(right=True, top=True)
        # plt.subplots_adjust(hspace=.5)
        plt.tight_layout()
        sns_plot.savefig("%s.pdf" % i_file.split("_")[0], dpi=400, bbox_inches='tight')


#################
def transform_data(in_file, hue_type):
    df = pd.read_table(in_file, sep="\t").stack()
    df_high, df_low = pd.DataFrame(df[:, "high_methy_cv"]), pd.DataFrame(df[:, "low_methy_cv"])
    df_high["term"], df_low["term"] = "HMC", "LMC"
    df_high.columns, df_low.columns = ["score", "term"], ["score", "term"]
    df_com = pd.DataFrame(pd.concat([df_low, df_high]))
    df_com["score"] = df_com["score"].astype(float)
    df_com["hue_type"] = hue_type
    print(df_com.head())
    return df_com


##################################################################################
def figure_2_c():   # #00BFFF #FF8C00
    os.chdir("/data3/xs/tissue_m6a/2018.1/fig2/fig2_315/1gene_filter_fpkm0")
    df = pd.read_table("Kmeans.txt", sep="\t").stack()
    df_high, df_low = pd.DataFrame(df[:, "high_methy_cv"]), pd.DataFrame(df[:, "low_methy_cv"])
    df_high["term"], df_low["term"] = "HMC", "LMC"
    df_high.columns, df_low.columns = ["score", "term"], ["score", "term"]
    # print(df_high.head())
    df_com = pd.DataFrame(pd.concat([df_low, df_high]))
    df_com["score"] = df_com["score"].astype(float)
    print(df_com.head())
    df_com["Category"] = "LMC    HMC"
    sns.set_style("white")
    f, axes = plt.subplots(figsize=(3, 3))
    flatui = ["#00BFFF", "#FF8C00"]
    sns_plot = sns.violinplot(x=df_com["Category"], y=df_com["score"], hue=df_com["term"], split=True, palette=sns.color_palette(flatui), width=.4, inner="quart").get_figure()    #  scale="area", width=1.0 palette=sns.color_palette(flatui),
    axes.tick_params(direction='out', labelsize=20)
    axes.tick_params(axis="x", labelsize=0)
    axes.set_xlabel("")
    axes.set_ylabel("")
    l = axes.legend()
    l.set_title("")
    # axes.set_ylabel("Relative expression\ndynamics of gene", fontsize=24, weight='bold', labelpad=12)
    # 去掉多余的线
    sns.despine(right=True, top=True)
    # plt.subplots_adjust(hspace=.5)
    plt.tight_layout()
    sns_plot.savefig("/data5/galaxy/project/plot_figure/Kmeans.pdf", dpi=400, bbox_inches='tight')


##################################################################################
def figure_2_d():   # #00BFFF #FF8C00
    os.chdir("/data3/xs/tissue_m6a/2018.1/fig2/fig2_315/1gene_filter_fpkm0")
    df_m = transform_data("Kmeans_mRNA.txt", "mRNA")
    df_l = transform_data("Kmeans_lincRNA.txt", "lincRNA")
    df_com = pd.DataFrame(pd.concat([df_m, df_l]))
    sns.set_style("white")
    f, axes = plt.subplots(figsize=(5, 3))
    flatui = ["#00BFFF", "#FF8C00"]
    sns_plot = sns.violinplot(x=df_com["hue_type"], y=df_com["score"], hue=df_com["term"], palette=sns.color_palette(flatui), width=.4, legend=False).get_figure()    #  scale="area", width=1.0 palette=sns.color_palette(flatui),
    axes.tick_params(direction='out', labelsize=20)
    axes.tick_params(axis="x", labelsize=0)
    axes.set_xlabel("")
    axes.set_ylabel("")
    # l = axes.legend()
    plt.legend(loc="best")
    # l.set_title("")
    # axes.set_ylabel("Relative expression\ndynamics of gene", fontsize=24, weight='bold', labelpad=12)
    # 去掉多余的线
    sns.despine(right=True, top=True)
    # plt.subplots_adjust(hspace=.5)
    plt.tight_layout()
    sns_plot.savefig("/data5/galaxy/project/plot_figure/Kmeans_2d.pdf", dpi=400, bbox_inches='tight')


# figure_2_d()


##################################################################################
def figure_3_d():   # #00BFFF #FF8C00
    os.chdir("/data3/xs/mettl3chip_hela")
    df_a = pd.read_table("mettl3_promoter_cpg", sep="\t", header=None, names=["score"]).drop_duplicates().dropna(how="any")
    df_a["term"] = "+"
    df_b = pd.read_table("no_mettl3_promoter_cpg", sep="\t", header=None, names=["score"]).drop_duplicates().dropna(how="any")
    df_b["term"] = "-"
    df_com = pd.DataFrame(pd.concat([df_b, df_a]))
    sns.set_style("white")
    f, axes = plt.subplots(figsize=(4, 4))
    flatui = ["#00BFFF", "#FF8C00"]
    sns_plot = sns.violinplot(x=df_com["term"], y=df_com["score"], palette=sns.color_palette(flatui), width=.4, legend=False).get_figure()    #  scale="area", width=1.0 palette=sns.color_palette(flatui),
    axes.tick_params(direction='out', labelsize=20)
    axes.tick_params(axis="x", labelsize=20)
    axes.set_xlabel("")
    axes.set_ylabel("")
    # l = axes.legend()
    plt.legend(loc="best")
    # l.set_title("")
    # axes.set_ylabel("Relative expression\ndynamics of gene", fontsize=24, weight='bold', labelpad=12)
    # 去掉多余的线
    sns.despine(right=True, top=True)
    # plt.subplots_adjust(hspace=.5)
    plt.tight_layout()
    sns_plot.savefig("/data5/galaxy/project/plot_figure/Kmeans_mettl3.pdf", dpi=400, bbox_inches='tight')


#
# figure_3_d()
figure_2_a_b()