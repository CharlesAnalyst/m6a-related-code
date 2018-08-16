#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob


m6a_dir = "/data5/galaxy/project/data/total_m6a_peak"
peak_fa_list = glob.glob("%s/*.fa" % m6a_dir)
result_dir = m6a_dir


def fasta_get_markov(i_fasta):
    i_markov = os.path.join(result_dir, "%s.txt" % os.path.basename(i_fasta).split(".fa")[0])
    os.system("fasta-get-markov %s %s" % (i_fasta, i_markov))


def main_method():
    for peak_fa in peak_fa_list:
        fasta_get_markov(peak_fa)

    os.system("cat ")


if __name__ == '__main__':
    main_method()
