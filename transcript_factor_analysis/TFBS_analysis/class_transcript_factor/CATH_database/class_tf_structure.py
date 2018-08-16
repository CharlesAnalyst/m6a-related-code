#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


data_file = "/data5/galaxy/project/tf_analysis/class_tf/CATH_db/cath-domain-description-file.txt"
result_file = "/data5/galaxy/project/tf_analysis/class_tf/CATH_db/cluster_result_structure.txt"

pdb_dict = {}

with open(data_file, 'r') as f:
    contents = f.readlines()
    i = 0
    while i < len(contents):
        if contents[i].startswith("DOMAIN"):
            domain = contents[i].split("DOMAIN")[1].strip()
            sub_dict = {}
            pdb_dict[domain] = sub_dict
            i += 1
            while i < len(contents) and not contents[i].startswith("DOMAIN"):
                if contents[i].startswith("CLASS"):
                    pdb_dict[domain]["class"] = contents[i].split("CLASS")[1].strip()
                elif contents[i].startswith("ARCH"):
                    pdb_dict[domain]["arch"] = contents[i].split("ARCH")[1].strip()
                elif contents[i].startswith("TOPOL"):
                    pdb_dict[domain]["topol"] = contents[i].split("TOPOL")[1].strip()
                elif contents[i].startswith("HOMOL"):
                    pdb_dict[domain]["homol"] = contents[i].split("HOMOL")[1].strip()
                i += 1
        else:
            i += 1

print(pdb_dict)
