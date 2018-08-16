#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
from multiprocessing import Pool


############################
species = "mouse"
TOTAL_THREAD_NUM = 30
#####################
gene_anno, bam_list, result_dir, merge_dir = "", [], "", ""
if species == "mouse":
    gene_anno = "/data/database/GRCm38/GENCODE/gencode.vM15.annotation.gff3"
    bam_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/unique_bam"
    bam_list = [x for x in glob.glob("%s%s*_sorted_unique.bam" % (bam_dir, os.sep)) if "_IP_" not in x]
    result_dir = "/data5/galaxy/project/DNMT1_KO/RNA_Seq/input/expression/stringtie"
    merge_dir = os.path.join(result_dir, "02.merge_assembly")
elif species == "human":
    gene_anno = "/data4/database/GRCh38/GENCODE/gencode.v27.annotation.gff3"
    bam_dir = "/data5/galaxy/project/unique_bam/input"
    bam_list = glob.glob("%s%s*_sorted_unique.bam" % (bam_dir, os.sep))
    result_dir = "/data5/galaxy/project/expression/stringtie/input"
    merge_dir = os.path.join(result_dir, "02.merge_assembly")
# ###########
for i_result_dir in [result_dir, merge_dir]:
    if not os.path.exists(i_result_dir):
        os.makedirs(i_result_dir)
each_thread = TOTAL_THREAD_NUM / len(bam_list)
#######################################################################################


def get_prefix_from_bam(bam):
    # contain_words = "_sorted_unique.bam"
    prefix = os.path.basename(bam).split("_")[0]
    return prefix


def first_assemble(bam):
    prefix = get_prefix_from_bam(bam)
    prime_dir = os.path.join(result_dir, prefix, "01.prime_assembly")
    if not os.path.exists(prime_dir):
        os.makedirs(prime_dir)
    out_transcripts = os.path.join(prime_dir, "assembled_transcripts.gtf")
    out_gene_abund = os.path.join(prime_dir, "gene_abund.tab")
    #  -m 50
    os.system("stringtie %s -G %s -p %d --rf -o %s -A %s -B -e" % (bam, gene_anno, each_thread, out_transcripts, out_gene_abund))


def merge_assemble():
    # Generate the primary gtf list file
    gtf_list = [os.path.abspath(gtf) for gtf in glob.glob("%s/*/01.prime_assembly/*.gtf" % result_dir)]
    merged_gtf = os.path.join(merge_dir, "merged_transcripts.gtf")
    # set minimum FPKM value  -F %f
    os.system("stringtie --merge -o %s %s" % (merged_gtf, " ".join(gtf_list)))


def second_assemble(bam):
    prefix = get_prefix_from_bam(bam)
    final_dir = os.path.join(result_dir, prefix, "03.final_assembly")
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
    new_gtf = os.path.join(merge_dir, "merged_transcripts.gtf")
    out_transcripts = os.path.join(final_dir, "assembled_transcripts.gtf")
    out_gene_abund = os.path.join(final_dir, "gene_abund.tab")
    #  -m 50
    os.system("stringtie %s -G %s -p %d --rf -o %s -A %s -B -e" % (bam, new_gtf, each_thread, out_transcripts, out_gene_abund))

"""
def filter_result(final_dir):
    in_transcripts = os.path.join(final_dir, "t_data.ctab")
    out_transcripts = os.path.join(final_dir, "t_data_fpkm-%s.ctab" % min_fpkm)
    data = pd.read_table(in_transcripts, sep="\t")
    sub_data = data[data["FPKM"] > float(min_fpkm)]
    sub_data.to_csv(out_transcripts, sep="\t", index=False)
    # ##
    in_gtf = os.path.join(final_dir, "assembled_transcripts.gtf")
    out_gtf = os.path.join(final_dir, "assembled_transcripts_fpkm-%s.gtf" % min_fpkm)
    with open(out_gtf, 'w') as fw:
        with open(in_gtf, 'r') as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    data_type = line.split("\t")[2]
                    if data_type == "transcript":
                        fpkm = float(line.split("\t")[8].split('FPKM "')[1].split('";')[0])
                        if fpkm > float(min_fpkm):
                            fw.write(line)
"""

if __name__ == "__main__":
    #
    pool = Pool(processes=len(bam_list))
    for i_bam in bam_list:
        pool.apply_async(first_assemble, (i_bam, ))
    pool.close()
    pool.join()
    #
    merge_assemble()
    #
    pool_2 = Pool(processes=len(bam_list))
    for i_bam in bam_list:
        pool_2.apply_async(second_assemble, (i_bam, ))
    pool_2.close()
    pool_2.join()
