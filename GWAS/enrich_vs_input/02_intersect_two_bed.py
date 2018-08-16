#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
import os
import glob
import datetime
from multiprocessing import Pool


m6a_bed_dir = "/data5/galaxy/project/data/total_m6a_peak/autosomal_peak"
con_bed_dir = "/data5/galaxy/project/data/input_data/autosomal/100cycle_input"
snp_dir = "/data5/galaxy/project/GWAS_db/04_class_snp/snp_classification"
result_dir = "/data5/galaxy/project/GWAS_analysis/intersect_snp_class"   #
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def get_input_list():
    os.chdir(con_bed_dir)
    bed_list = glob.glob("control-*_*.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_tissue_name_from_input(input_file):
    tissue_name = os.path.basename(input_file).split("-")[1].split(".bed")[0].lower()
    return tissue_name


def get_m6a_list():
    os.chdir(m6a_bed_dir)
    bed_list = glob.glob("*.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_tissue_from_m6a(m6a_file):
    tissue_name = os.path.basename(m6a_file).split(".bed")[0].lower()
    return tissue_name


def get_snp_list():
    os.chdir(snp_dir)
    bed_list = glob.glob("gwas_snp_*.bed")
    # bed_list = glob.glob("*.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_snp_name_from_snp(snp_file):
    snp_name = os.path.basename(snp_file).split("_")[2].split(".bed")[0]
    # snp_name = os.path.basename(snp_file).split(".bed")[0]
    return snp_name


def intersect_bed_snp(peak_type, peak_bed, snp_bed):
    snp_name = get_snp_name_from_snp(snp_bed)
    intersect_result_file = ""
    if peak_type == "m6a":
        m6a_name = get_tissue_from_m6a(peak_bed)
        intersect_result_file = "%s_%s.bed" % (snp_name, m6a_name)
    elif peak_type == "control":
        con_name = get_tissue_name_from_input(peak_bed)
        intersect_result_file = "%s_%s.bed" % (snp_name, con_name)
    os.system("bedtools intersect -a {snp_bed} -b {peak_bed} -wa > {intersect_result_file}"
              .format(snp_bed=snp_bed, peak_bed=peak_bed, intersect_result_file=intersect_result_file))
    print("%s\t%s done!" % (peak_type, snp_name))


def run_all():
    m6a_list = get_m6a_list()
    input_list = get_input_list()
    snp_list = get_snp_list()
    print(len(m6a_list), len(input_list), len(snp_list))
    #
    os.chdir(m6a_bed_dir)
    condition_result_dir = "%s%scondition" % (result_dir, os.sep)
    if not os.path.exists(condition_result_dir):
        os.makedirs(condition_result_dir)
    os.chdir(condition_result_dir)
    for m6a in m6a_list:
        for snp in snp_list:
            intersect_bed_snp("m6a", m6a, snp)
    #
    os.chdir(con_bed_dir)
    control_result_dir = "%s%scontrol" % (result_dir, os.sep)
    if not os.path.exists(control_result_dir):
        os.makedirs(control_result_dir)
    os.chdir(control_result_dir)
    pool = Pool(processes=20)
    for con in input_list:
        for snp in snp_list:
            pool.apply_async(intersect_bed_snp, ("control", con, snp))
    pool.close()
    pool.join()


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_all()
    print(start_time)
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
