'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-05 11:02:33
LastEditors: Ne0tea
LastEditTime: 2023-11-05 13:24:20
'''
import difflib
import sys

gene_file=sys.argv[1]
anno_file=sys.argv[2]
gene_list=[]
with open(gene_file,'r') as gf:
    for i in gf:
        gene_list.append(i.strip())

with open(anno_file,'r') as af:
    for i in af:
        line=i.strip().split()
        if line[3]=='gene':
            range=[int(line[1]),int(line[2])]
            cutoff=int(line[1])
            start=0
            end=int(line[2])-cutoff
        else:
            start=int(line[1])-cutoff
            end=int(line[2])-cutoff
        try:
            gene_name=difflib.get_close_matches(line[0].split('=')[1],gene_list,1, cutoff=0.4)[0]
        except IndexError:
            print('Suspect',line[0].split('=')[1],'in',gene_name)
        print(gene_name,start,end,line[3],line[4],sep='\t')
