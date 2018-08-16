#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

from scipy.stats import norm
import pandas as pd
import math


#########################################
# cdf: Cumulative Distribution Function;
# sf: Survival Function (1-CDF)
###################################################################################
in_file = "/data5/galaxy/project/tf_analysis/bg_distribution"
input_values = [0.1, 0.6, 1.2, 1.5]
result_file = "/data5/galaxy/project/tf_analysis/bg_distribution/pvalues_result.txt"
####################################################################################


def transform_to_zscore(value, mean, std, n):
    z_score = (mean - value) / (std / math.sqrt(n))
    return z_score


def calculate_pvalue_from_zscore(z_score):
    p_values = norm.sf(abs(z_score)) * 2  # two sided
    return p_values


def get_p_value():
    with open(result_file, 'w') as fw:
        fw.write("input_value\tp_value\n")
        df = pd.read_table(in_file, sep="\t", header=None)
        mean, std, n = df.mean(), df.std(), len(df)
        for value in input_values:
            z = transform_to_zscore(value, mean, std, n)
            pvalue = calculate_pvalue_from_zscore(z)
            fw.write("%f\t%f\n" % (value, pvalue))
