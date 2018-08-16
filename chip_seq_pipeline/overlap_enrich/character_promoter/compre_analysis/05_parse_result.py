#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import math
from scipy import stats
import numpy as np
import seaborn as sns
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

##################################################################
os.chdir("/data5/galaxy/project/CpG_m6a_motif/fasta_seq")
high_file, low_file = "high.count", "low.count"
###########################################################


def main():
    high_counts, low_counts = get_count_list(high_file), get_count_list(low_file)
    print(high_counts)
    print(low_counts)
    h_mean, low_mean = np.mean(high_counts), np.mean(low_counts)
    print(h_mean, low_mean)
    print(max(high_counts), max(low_counts))
    print(min(high_counts), min(low_counts))
    sns.boxplot(data=high_counts)
    sns.boxplot(data=low_counts)
    plt.show()
    # plot_histogram(high_counts)
    # plot_histogram(low_counts)
    statistic, pvalue = statistic_significance(high_counts, low_counts)
    print(statistic, pvalue)


def get_count_list(i_file):
    result_list = []
    with open(i_file, 'r') as f:
        for line in f.readlines():
            seq = line.split("{")[-1].split("}}}")[0]
            number = len(seq.split(","))
            result_list.append(number)
    return result_list


def plot_histogram(data_list):
    num_bins, mean, std = 20, np.mean(data_list), np.std(data_list)
    n, bins, patches = plt.hist(data_list, num_bins, normed=1, facecolor="blue", alpha=0.5)
    y = mlab.normpdf(bins, mean, std)
    plt.plot(bins, y, "r--")
    plt.xlabel("Number of RRACH occurence")
    plt.ylabel("Probability")
    plt.title(r"Histogram of IQ: $mean=%s$, $std=%s$" % (str(mean), str(std)))
    plt.subplots_adjust(left=0.15)
    plt.show()


def statistic_significance(m6a_scores, free_scores):
    statistic, pvalue = stats.mannwhitneyu(m6a_scores, free_scores)
    return statistic, pvalue


if __name__ == '__main__':
    main()