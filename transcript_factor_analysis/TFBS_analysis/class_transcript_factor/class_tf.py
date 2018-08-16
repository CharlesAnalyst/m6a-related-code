#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
import itertools
from multiprocessing import Pool


pdb_dir = "/data5/galaxy/project/tf_analysis/class_tf/PDB_data"
pdb_list = glob.glob("%s/*.pdb" % pdb_dir)
result_dir = "/data5/galaxy/project/tf_analysis/class_tf/ssap_result"
align_dir = "%s/align_result" % result_dir
if not os.path.exists(align_dir):
    os.makedirs(align_dir)
result = os.path.join(result_dir, "ssap_result.txt")


def cath_ssap_run(pdb_prefix_list):
    pdb_1, pdb_2 = pdb_prefix_list[0], pdb_prefix_list[1]
    pdb_1, pdb_2 = os.path.basename(pdb_1), os.path.basename(pdb_2)
    os.system("cath-ssap.ubuntu14.04 --aligndir %s --pdb-path %s %s %s >> %s" % (align_dir, pdb_dir, pdb_1, pdb_2, result))


if __name__ == '__main__':
    combination_list = itertools.combinations(pdb_list, 2)
    pool = Pool()
    pool.map_async(cath_ssap_run, combination_list)
    pool.close()
    pool.join()
