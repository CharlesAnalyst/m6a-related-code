#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import random
import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#####################################################################################
data_dir = "/data5/galaxy/project/promoter_TF_enrich/comprehensive_enrich/total_gene"
bin_num = 15
#################################################################


def main():
    tissue_dict = get_tissue_dict()
    for tissue, files in tissue_dict.items():
        print(tissue)
        get_data_distribution(files)


def get_tissue_dict():
    tissue_dict = {}
    file_list = glob.glob("%s/*.txt" % data_dir)  # overlap-number_1-kidney.txt
    for i_file in file_list:
        tissue = os.path.basename(i_file).split("-")[-1].split(".txt")[0].lower()
        tissue_dict[tissue] = tissue_dict.get(tissue, []) + [i_file]
    return tissue_dict


def get_data_distribution(tissue_files):
    real_group_list = get_real_bin_group(tissue_files)
    plot_values, histogram_range = stat_each_bin_sum(real_group_list)
    plot_histogram(plot_values, histogram_range)


def get_real_bin_group(tissue_files):
    real_group_list = []
    for i_file in tissue_files:
        data_list = pd.read_table(i_file, sep="\t", header=None, names=["gene", "number"]).loc[:, "number"].tolist()
        real_group_list.append(data_list)
    return real_group_list


def stat_each_bin_sum(real_group_list):
    combined_list, total_list, bin_edge_list, plot_values = [], [], [], []
    for i_list in real_group_list:
        combined_list += i_list
    histogram_range = (min(combined_list), max(combined_list))
    for i_list in real_group_list:
        num_list, bin_edge_list = stat_each_bin_number(i_list, histogram_range)
        total_list.append(num_list)
    df = pd.DataFrame(np.array(total_list))
    # print(df.head())
    sum_num_list = [df[col].sum() for col in df]
    for i in range(len(bin_edge_list)):
        random_value_list = random_generate_value(bin_edge_list[i][0], bin_edge_list[i][1], sum_num_list[i])
        plot_values += random_value_list
    return plot_values, histogram_range


def stat_each_bin_number(value_list, histogram_range):
    num_list, bins, patches = plt.hist(value_list, bin_num, range=histogram_range)
    bin_edge_list = []
    for i in range(len(bins)-1):
        if i == (len(bins) - 2):
            bin_edge = (bins[i], bins[i+1])
        else:
            bin_edge = (bins[i], (bins[i+1] - 0.00001))
        bin_edge_list.append(bin_edge)
    return num_list, bin_edge_list


def random_generate_value(lower, upper, value_num):
    random_value_list = []
    for i in range(int(value_num)):
        random_value_list.append(random.uniform(lower, upper))  # [lower, upper]
    return random_value_list


def plot_histogram(data_list, range_tuple):
    num_bins, mean, std = 20, np.mean(data_list), np.std(data_list)
    n, bins, patches = plt.hist(data_list, num_bins, range=range_tuple, normed=1, facecolor="blue", alpha=0.5)
    y = mlab.normpdf(bins, mean, std)
    plt.plot(bins, y, "r--")
    plt.xlabel("Number of overlapped transcript factor")
    plt.ylabel("Probability")
    plt.title(r"Histogram of IQ: $mean=%s$, $std=%s$" % (str(mean), str(std)))
    plt.subplots_adjust(left=0.15)
    plt.show()


# data_list = [1, 2, 2, 3, 4, 4, 4, 6, 6, 6]
# n, bins, patches = plt.hist(data_list, 5, range=(1, 10))
# print(n[0])
if __name__ == '__main__':
    main()