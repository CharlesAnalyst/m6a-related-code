#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os


ssap_file = "/data5/galaxy/project/tf_analysis/class_tf/ssap_result/ssap_result.txt"
result_dir = "/data5/galaxy/project/tf_analysis/class_tf/ssap_result"
result_file = os.path.join(result_dir, "correlation_matrix.txt")


def read_file():
    with open(ssap_file, 'r') as f:
        contents = f.readlines()
    return contents


def get_total_protein(contents):
    pdb_1_list, pdb_2_list = [], []
    for line in contents:
        pdb_1_list.append(line.split()[0]), pdb_2_list.append(line.split()[1])
    total_protein = list(set(pdb_1_list + pdb_2_list))
    return total_protein


def get_score_dict(contents):
    pdb_dict = {}
    for line in contents:
        pdb_1, pdb_2, ssap_score = line.split()[0], line.split()[1], line.split()[2]
        pdb_dict["%s-%s" % (pdb_1, pdb_1)] = 1.0
        pdb_dict["%s-%s" % (pdb_2, pdb_2)] = 1.0
        pdb_dict["%s-%s" % (pdb_1, pdb_2)] = ssap_score
        pdb_dict["%s-%s" % (pdb_2, pdb_1)] = ssap_score
    return pdb_dict


def main_method():
    contents = read_file()
    total_protein = get_total_protein(contents)
    score_dict = get_score_dict(contents)
    title = "\t".join(total_protein)
    with open(result_file, 'w') as fw:
        fw.writelines("Protein\t%s\n" % title)
        for i in range(len(total_protein)):
            index_name, line_score_list = total_protein[i], []
            for j in range(len(total_protein)):
                col_name = total_protein[j]
                try:
                    line_score_list.append(str(score_dict["%s-%s" % (index_name, col_name)]))
                except KeyError:
                    print("can't find %s-%s" % (index_name, col_name))
            fw.writelines("%s\t%s\n" % (index_name, "\t".join(line_score_list)))


if __name__ == '__main__':
    main_method()
