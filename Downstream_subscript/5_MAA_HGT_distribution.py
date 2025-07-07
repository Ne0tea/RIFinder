'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-10-24 20:51:52
LastEditors: Ne0tea
LastEditTime: 2023-10-28 16:53:23
'''

import numpy as np
from pandas import Series,DataFrame
from collections import Counter

def threshold_cluster(Data_set,threshold):
    #统一格式化数据为一维数组
    stand_array=np.asarray(Data_set).ravel('C')
    stand_Data=Series(stand_array)
    index_list,class_k=[],[]
    while stand_Data.any():
        if len(stand_Data)==1:
            index_list.append(list(stand_Data.index))
            class_k.append(list(stand_Data))
            stand_Data=stand_Data.drop(stand_Data.index)
        else:
            class_data_index=stand_Data.index[0]
            class_data=stand_Data[class_data_index]
            stand_Data=stand_Data.drop(class_data_index)
            if (abs(stand_Data-class_data)<=threshold).any():
                args_data=stand_Data[abs(stand_Data-class_data)<=threshold]
                stand_Data=stand_Data.drop(args_data.index)
                index_list.append([class_data_index]+list(args_data.index))
                class_k.append([class_data]+list(args_data))
            else:
                index_list.append([class_data_index])
                class_k.append([class_data])
    return index_list,class_k

gene_len={}
gene_ord={}
bed_file=r'E:\Bio_analysis\Grass_horizatal_transversion\subg_gene.order'
with open (bed_file,'r') as bedf:
    for i in bedf.readlines():
        c_line=i.strip().split()
        gene_id=c_line[0]
        contig=c_line[1]
        order=c_line[2]
        sp=gene_id.split('_')[0]
        # lennn=int(gene_e)-int(gene_s)
        # gene_len[gene_id]=lennn
        if sp not in gene_ord:
            gene_ord[sp]={gene_id:order}
        else:
            gene_ord[sp][gene_id]=order
# print(gene_ord)
HGT_gene_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\transfer_summary.txt"

hgt_gene=[]
og_dic={}
with open(HGT_gene_file,'r') as hgt_file:
    for i in hgt_file:
        c_gene=i.strip().split(' ')[1]
        og_dic[c_gene]=i.strip().split(' ')[0]
        hgt_gene.append(c_gene)
subg_list=set([x.split('_')[0] for x in hgt_gene])
block=0
block_dic={}
block_len_set=[]
for i in subg_list:
    ord_list=[int(gene_ord[i][x]) for x in hgt_gene if x.split('_')[0]==i]
    gene_list=[x for x in hgt_gene if x.split('_')[0]==i]
    index_list,class_k=threshold_cluster(ord_list,10)
    # print(class_k)
    for y in index_list:
        block+=1
        block_gene=[gene_list[x] for x in y]
        # if 'TaesA_TraesCS2A01G051700' in block_gene:
        #     print(block_gene)
        block_dic['block'+str(block)]=block_gene
        block_len_set.append(len(block_gene))
result = Counter(block_len_set)
print(result)

gff_file=r'E:\Bio_analysis\Grass_horizatal_transversion\ful_gff.bed'
gene_start={}
gene_end={}
with open (gff_file,'r') as gfff:
    for i in gfff.readlines():
        c_line=i.strip().split()
        gene_id=c_line[1]
        contig=c_line[0]
        start=min(int(c_line[2]),int(c_line[3]))
        end=max(int(c_line[2]),int(c_line[3]))
        sp=gene_id.split('_')[0]
        gene_start[gene_id]=start
        gene_end[gene_id]=end

out_file=open('HGT_dis.out','w')
large_block=[block_dic[x] for x in block_dic if len(block_dic[x]) >= 2]
for i in large_block:
    block_start=min([gene_start[x] for x in i])
    block_end=max([gene_end[x] for x in i])
    block_og=[og_dic[x] for x in i]
    if len(set(block_og)) >1:
        line="|".join(i)+"\t"+str(len(i))+"\t"+"|".join(block_og)+"\t"+str(int(block_end)-int(block_start))+"\t"+str(block_start)+"\t"+str(block_end)
        out_file.write(line+'\n')