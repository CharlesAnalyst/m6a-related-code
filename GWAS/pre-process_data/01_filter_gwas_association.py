#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy
# dbSNP Build 149, Genome Assembly GRCh38.p10, NCBI
# exclude SNPs mapping outside the main chromosome
# exclude SNPs without coordinates in GRCh38.p10 human genome assembly
# exclude SNPs without a dbSNP ID
# exclude records which were a combination of multiple SNPs associated with a disease or trait.（多个snp在同一行）
# get unique SNP-disease/trait combinations;


gwas_association_file = "/data5/galaxy/project/GWAS_db/gwas-catalog-associations_ontology-annotated.tsv"
result_file = "/data5/galaxy/project/GWAS_db/filtered_association.tsv"
main_chr_list = ["x", "y"]
for i in range(23):
    main_chr_list.append(str(i))


def read_file():
    with open(gwas_association_file, 'r', encoding='utf-8') as f:
        title = f.readline().strip()
        contents = f.readlines()
    return title, contents


def filter_file(contents):
    result_list = []
    for line in contents:
        info = line.strip().split("\t")
        chr_id, chr_pos, rs_id, mapped_trait = info[11], info[12], info[21], info[34]
        # print(chr_id, chr_pos, rs_id, mapped_trait)
        if chr_id.lower() in main_chr_list:
            if chr_pos != "":
                if rs_id.startswith("rs") and ";" not in rs_id:
                    if mapped_trait != "":
                        result_list.append(line.strip())
    return result_list


def split_mul_trait(result_list):
    tmp_list = []
    for line in result_list:
        info = line.strip().split("\t")
        mapped_trait = info[34]
        prefix = line.strip().split(mapped_trait)[0]
        postfix = line.strip().split(mapped_trait)[1]
        if "," in mapped_trait:
            snp_list = mapped_trait.split(", ")
            for snp in snp_list:
                tmp_list.append("%s%s%s" % (prefix, snp, postfix))
        else:
            tmp_list.append(line.strip())
    return tmp_list


def write_to_file(title, result_list):
    with open(result_file, "w", encoding='utf-8') as fw:
        fw.write("%s\r" % title.strip())
        for line in result_list:
            fw.write("%s\r" % line.strip())


if __name__ == "__main__":
    file_title, file_contents = read_file()
    file_result_list = filter_file(file_contents)
    final_result_list = split_mul_trait(file_result_list)
    write_to_file(file_title, final_result_list)
