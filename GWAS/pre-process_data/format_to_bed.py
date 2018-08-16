#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

# Generate GRCH38(1-base) version
csv_file = "I:\\PROJECT\\m6A_R-loop_snp\\SNP\\GWAS\\replicated_internal.txt"
out_file = "I:\\PROJECT\\m6A_R-loop_snp\\SNP\\GWAS\\replicated_internal.bed"


def read_tsv():
    with open(csv_file, 'r') as f:
        result_list = []
        title = f.readline().strip()
        for line in f.readlines():
            info = line.strip().split("\t")
            chr_id, chr_pos, rs_id = ("chr%s" % info[11]), info[12], info[21]
            chr_start = str(int(chr_pos) - 1)
            chr_end = chr_pos
            result_list.append("%s\t%s\t%s\t%s" % (chr_id, chr_start, chr_end, rs_id))
        return result_list


def write_to_file(result_list):
    result_list = list(set(result_list))
    result_list.sort(key=lambda d: (d.split("\t")[0], int(d.split("\t")[1])))
    with open(out_file, 'w') as fw:
        for line in result_list:
            fw.write("%s\n" % line.strip())


if __name__ == "__main__":
    tsv_result_list = read_tsv()
    write_to_file(tsv_result_list)