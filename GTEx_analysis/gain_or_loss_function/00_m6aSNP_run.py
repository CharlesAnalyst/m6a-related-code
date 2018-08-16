#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob
from multiprocessing import Pool


######################################################################################################################
m6aSNP_jar = "/data/software/m6ASNP/m6ASNP.jar"
eQTL_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/eQTL_snp_db/GTEx_Analysis_v7_eQTL/01_1_format_vcf_hg19"
ucsc_annotation = "/data/software/m6ASNP/knownGeneAnnotation_hg19"
genome_sequence = "/data/software/m6ASNP/hg19.2bit"
result_dir = "/data5/galaxy/project/snp_analysis/GTEx_analysis/gain_or_loss/01_predict_result_hg19"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
###################################


def predict_single_file(eQTL_vcf):
    result_txt = os.path.join(result_dir, os.path.basename(eQTL_vcf).split(".txt")[0])
    os.system("java -jar %s -predict -i %s -it tab -sp Human -a %s -g %s -o %s" %
              (m6aSNP_jar, eQTL_vcf, ucsc_annotation, genome_sequence, result_txt))
    # java -jar /data/software/m6ASNP/m6ASNP.jar -predict -i liver.txt -it tab -sp Human -a /data/software/m6ASNP/knownGeneAnnotation_hg19 -g /data/software/m6ASNP/hg19.2bit -o test_result


if __name__ == '__main__':
    eQTL_vcf_list = glob.glob("%s/*.txt" % eQTL_dir)
    pool = Pool()
    for eQTL_vcf in eQTL_vcf_list:
        pool.apply_async(predict_single_file, (eQTL_vcf, ))
    pool.close()
    pool.join()

"""
# gwas snp annotate
gwas_bed = "/data5/galaxy/project/GWAS_db/snp_classification/merge/total_gwas_hg38.bed"
ucsc_annotation = "/data/software/m6ASNP/knownGeneAnnotation_hg38"
genome_sequence = "/data/software/m6ASNP/hg38.2bit"
os.system("java -jar %s -predict -i %s -it tab -sp Human -a %s -g %s -o %s" %
            (m6aSNP_jar, eQTL_vcf, ucsc_annotation, genome_sequence, result_txt))
"""