#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


def calculate_pvalue(query_value, db_values):
    db_values = [float(x) for x in db_values]
    right_values = [x for x in db_values if x > float(query_value)]
    pvalue = len(right_values) / (len(db_values) + 1)
    print(pvalue)
    return pvalue