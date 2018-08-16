#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
import os
import glob
import datetime
import pandas as pd
from scipy import stats
from multiprocessing import Pool


m6a_bed = "/data5/galaxy/project/data/total_m6a_peak/universal/universal.bed"
control_bed = "/data5/galaxy/project/data/input_data/autosomal/universal/input_bed/universal_discrete.bed"
gwas_snp_bed = "/data5/galaxy/project/GWAS_db/03_replicated_GRCh38/autosomol_gwas_snp.bed"
cycle_dir = "/data5/galaxy/project/data/input_data/autosomal/universal/100cycle_input"


def intersect_control_snp(i_bed):
    print(i_bed)
    control_popen = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (gwas_snp_bed, i_bed))
    control_num = control_popen.read().split()[0]
    return control_num


def intersect_bed_snp():
    m6a_popen = os.popen("bedtools intersect -a %s -b %s -wa | wc -l" % (gwas_snp_bed, m6a_bed))
    control_list = glob.glob("%s/*.bed" % cycle_dir)
    result_list = []
    pool = Pool()
    for i_bed in control_list:
        i_result = pool.apply_async(intersect_control_snp, (i_bed, ))
        result_list.append(i_result)
    pool.close()
    pool.join()
    m6a_num = int(m6a_popen.read().split()[0])
    control_num_list = [int(i_result.get()) for i_result in result_list]
    return m6a_num, control_num_list


def stats_fisher_exact(m6a_snp, con_snp):
    m6a_all, con_all = statics_bed_sum(m6a_bed), statics_bed_sum(control_bed)
    m6a_unsnp, con_unsnp = (m6a_all - m6a_snp), (con_all - con_snp)
    #        m6a   control
    # snp    a       b
    # unsnp  c       d
    if m6a_snp == con_snp == 0 or m6a_unsnp == con_unsnp == 0 or m6a_snp == m6a_unsnp == 0 or con_snp == con_unsnp == 0:
        oddsratio, pvalue = 0, 999
    else:  # con_snp == 0 or m6a_unsnp == 0:  # Haldane-Anscombe correction#
        oddsratio, pvalue = stats.fisher_exact([[m6a_snp + 0.5, con_snp + 0.5], [m6a_unsnp + 0.5, con_unsnp + 0.5]])
    print(oddsratio, pvalue)
    if float(pvalue) > 0.05:
        oddsratio = 0
    return oddsratio


def statics_bed_sum(input_bed):
    df = pd.read_table(input_bed, sep="\t", header=None, comment="#", names=["chr", "start", "end"])
    df.loc[:, "length"] = df["end"] - df["start"]
    return df["length"].sum()


def run_100_fisher_test(m6a_num, control_num_list):
    pool, result_list = Pool(), []
    for control_num in control_num_list:
        oddsratio = pool.apply_async(stats_fisher_exact, (m6a_num, control_num))
        result_list.append(oddsratio)
    pool.close()
    pool.join()
    results = [float(x.get()) for x in result_list]
    mean_oddsratio = sum(results) / len(results)
    return mean_oddsratio


def main_method():
    m6a_num, control_num_list = intersect_bed_snp()
    mean_oddsratio = run_100_fisher_test(m6a_num, control_num_list)
    print("mean_oddsratio: %f" % mean_oddsratio)


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    main_method()
    print(start_time)
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
