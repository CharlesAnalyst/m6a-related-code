#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import datetime
import os


# map_file = "/data5/galaxy/project/GWAS_db/cluster_result-2018_1_19.txt"
map_file = "/data5/galaxy/project/GWAS_db/04_class_snp/cluster_by_map/cluster_result.txt"
filtered_file = "/data5/galaxy/project/GWAS_db/03_replicated_GRCh38/replicated_GWS_snp.txt"
result_dir = "/data5/galaxy/project/GWAS_db/04_class_snp/cluster_by_map/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def read_map():
    map_dict = {}
    with open(map_file, 'r') as f:
        contents = f.readlines()
        for line in contents:
            info = line.split("\t")
            parent_term, disease_trait, disease_num,  = info[0], info[1].lower(), info[2]
            if disease_trait not in map_dict:
                map_dict[disease_trait] = parent_term.strip()
    return map_dict


def read_disease_file():
    disease_list = []
    with open(filtered_file, 'r') as f:
        f.readline()
        for line in f.readlines():
            info = line.strip().split("\t")
            chr_id, chr_pos, rs_id, mapped_trait = ("chr%s" % info[11]), info[12], info[21], info[34].lower()
            chr_start, chr_end = str(int(chr_pos) - 1), chr_pos     # GWAS Catalog us
            disease_list.append("%s\t%s\t%s\t%s\t%s" % (chr_id, chr_start, chr_end, rs_id, mapped_trait))
    return disease_list


def query_relation_according_map():
    result_list = []
    disease_list = read_disease_file()
    map_dict = read_map()
    for line in disease_list:
        info = line.strip().split("\t")
        chr_id, chr_start, chr_end, rs_id, trait = info
        if trait in map_dict:
            parent_term = map_dict[trait]
        else:
            print(trait)
            continue
        new_line = "%s\t%s\t%s\t%s\t%s" % (parent_term, chr_id, chr_start, chr_end, rs_id)
        result_list.append(new_line)
    return result_list


def write_to_single_file(result_list, result_file):
    # 第一列是term
    result_list.sort(key=lambda d: (d.split("\t")[0], d.split("\t")[1].split("chr")[1], int(d.split("\t")[2])))
    unique_list = list(set(result_list))
    # count_list = [result_list.count(i) for i in unique_list]
    # combine_list = ["%s\t%d" % (unique_list[i], count_list[i]) for i in range(len(unique_list))]
    # combine_list.sort(key=lambda d: (d.split("\t")[0], int(d.split("\t")[2])), reverse=True)
    with open(result_file, "w", encoding='utf-8') as fw:
        for line in unique_list:
            fw.write("%s\r" % line.strip())
        # for i in range(len(combine_list)):
        #     fw.write("%s\r" % (combine_list[i].strip()))


def write_to_multi_file(result_list):
    os.chdir(result_dir)
    # print(list(set([line.split("\t")[0] for line in result_list])))
    parent_term_unique_list = list(set([line.split("\t")[0].replace(" ", "-") for line in result_list]))
    print(parent_term_unique_list)
    result_dict = {}
    for par_term in parent_term_unique_list:
        result_dict[par_term] = []
    for line in result_list:
        info = line.strip().split("\t")
        new_line = "\t".join(info[1:])
        parent_term = info[0].replace(" ", "-")
        try:
            result_dict[parent_term].append(new_line)
        except KeyError:
            print("%s doesn't exit in result dict" % parent_term)
    write_to_file(result_dict)


def write_to_file(result_dict):
    for term in result_dict:
        snp_list = result_dict[term]
        uniq_snp_list = list(set(snp_list))
        uniq_snp_list.sort(key=lambda d: (d.split("\t")[0].split("chr")[1], int(d.split("\t")[1])))
        outfile = "gwas_snp_%s.bed" % term
        with open(outfile, 'w') as fw:
            for line in uniq_snp_list:
                fw.write("%s\n" % line.strip())


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    all_result_list = query_relation_according_map()
    # write_to_single_file(all_result_list)
    write_to_multi_file(all_result_list)
    print(start_time)
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
