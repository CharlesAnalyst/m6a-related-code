#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import datetime
from multiprocessing import Pool
# http://www.bio-info-trainee.com/1211.html
# no markov file!

######################################################################################
ip = "/data5/galaxy/project/tf_pathway/m6a_snp_intersect/peak/all_tissues.bed"
# ip_list = glob.glob("%s/*.bed" % ip_dir)
reference_genome = "/data4/database/GRCh38/GENCODE/GRCh38.primary_assembly.genome.fa"
# motif_file = "/data4/database/tf_motif/motif_matrix.txt"
motif_file = "/data5/galaxy/project/TF_pathway_each_tissue/TF_db/motif_matrix.txt"
####################
result_dir = "/data5/galaxy/project/tf_pathway/fimo_result"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#######################################################################################


def identify_motif_from_fa(bed_file):
    fasta_file = "%s.fa" % bed_file.split(".bed")[0]
    os.system("bedtools getfasta -fi %s -bed %s -fo %s" % (reference_genome, bed_file, fasta_file))
    i_out_dir = os.path.join(result_dir, os.path.basename(fasta_file).split(".fa")[0])
    if not os.path.exists(i_out_dir):
        os.makedirs(i_out_dir)
    os.system("fimo --oc %s %s %s" % (i_out_dir, motif_file, fasta_file))


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # pool = Pool()
    # for ip in ip_list:
        # pool.apply_async(identify_motif_from_fa, (ip,))
    identify_motif_from_fa(ip)
    # pool.close()
    # pool.join()
    print(start_time)
    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(end_time)
