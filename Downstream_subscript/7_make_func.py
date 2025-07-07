'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-06 22:29:55
LastEditors: Ne0tea
LastEditTime: 2023-11-07 15:03:19
'''
import sys

hgt_file=sys.argv[1]
pfam_f=sys.argv[2]

gene_list=[]
with open(hgt_file,'r') as hf:
    for i in hf:
        gene_list.append(i.strip())

in_hgt=[]
with open(pfam_f,'r') as pf:
    for i in pf:
        traget_g=i.split()[0].replace('.','_')
        if traget_g in gene_list:
            print(i.strip())
            in_hgt.append(i.split()[0])
print(set(gene_list)-set(in_hgt))