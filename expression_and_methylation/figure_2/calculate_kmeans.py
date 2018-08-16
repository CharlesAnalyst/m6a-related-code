#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import numpy as np
from sklearn.cluster import KMeans


def k_means(k, value_list):
    value_list, high_index, median_index, low_index, expre_high, expre_median, expre_low = list(value_list), [], [], [], [], [], []
    # print(value_list)
    values = [[float(i)] for i in value_list if str(i) != "nan"]
    estimator = KMeans(n_clusters=k, random_state=0).fit(values)
    labels = estimator.labels_
    centers = estimator.cluster_centers_
    # sort labels, ascending
    center_list, label_list = [], []
    for i in range(len(centers)):
        center_list.append(centers[i][0])
        label_list.append(i)
    sorted_labels = [label_list[i] for i in np.argsort(center_list)]
    if k == 2:
        low_index = [i for i in range(len(labels)) if labels[i] == sorted_labels[0]]
        high_index = [i for i in range(len(labels)) if labels[i] == sorted_labels[1]]
        return high_index, low_index
    elif k == 3:
        low_index = [i for i in range(len(labels)) if labels[i] == sorted_labels[0]]
        median_index = [i for i in range(len(labels)) if labels[i] == sorted_labels[1]]
        high_index = [i for i in range(len(labels)) if labels[i] == sorted_labels[2]]
        return high_index, median_index, low_index

