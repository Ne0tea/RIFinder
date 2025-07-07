'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-01-01 15:25:00
LastEditors: Ne0tea
LastEditTime: 2025-01-01 16:06:49
'''
import sys

asm_file=sys.argv[1]
with open(asm_file,'r') as af:
	species_list = []
	for i in af:
		line = i.strip().split('\t')
		species_list.append(line[0])
species_set = set(species_list)

# 假设您的文件内容保存在一个文本文件中
with open(sys.argv[2], 'r') as f:
    lines = f.readlines()

missing_counts = {species: 0 for species in species_set}

for line in lines:
    species_in_line = {entry.split('_')[0] for entry in line.strip().split('\t')[4:]}

    for species in species_set:
        if species not in species_in_line:
            missing_counts[species] += 1

thresholds=[100,150,200,300,400,500]
results = {}
for threshold in thresholds:
    # 获取超过当前阈值的键
    keys_above_threshold = [key for key, value in missing_counts.items() if value > threshold]
    # 统计数量并保存
    results[threshold] = {
        "count": len(keys_above_threshold),
        "keys": keys_above_threshold,
    }

keep_sp = [ x for x in species_set if x not in results[150]['keys']]
# print(keep_sp)
for line in lines:
    species_in_line = {entry.split('_')[0] for entry in line.strip().split('\t')[4:]}
    species_in_line = set(species_in_line)
    if all([True if x in species_in_line else False for x in keep_sp]):
        for sp_gene in line.strip().split('\t')[4:]:
            cur_sp_gene = sp_gene.split(',')
            if cur_sp_gene[0].split('_')[0] in keep_sp:
                cur_gene_name = cur_sp_gene[0].split('_',1)[1]
                print(cur_gene_name)
    # else:
    #     print(species_in_line)