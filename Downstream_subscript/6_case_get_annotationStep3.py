'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-05 11:02:33
LastEditors: Ne0tea
LastEditTime: 2023-11-07 10:20:43
'''
from ete3 import Tree
import sys

gene_file=sys.argv[1]
# gene_file=r'E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\OG296_PP.list'
tree_file=sys.argv[2]
# tree_file=r'E:\Bio_analysis\Grass_horizatal_transversion\poaceae_sptree_rooted.nwk'
gene_list=[]
gene_dic={}
with open(gene_file,'r') as gf:
    for i in gf:
        gene_list.append(i.strip())

gene_tree=Tree(tree_file)
subtree_taxa = [x.split('_')[0] for x in gene_list]
for i in gene_list:
    if i.split('_')[0] in gene_dic:
        gene_dic[i.split('_')[0]].append(i)
    else:
        gene_dic[i.split('_')[0]]=[i]
gene_tree.prune(subtree_taxa)
# print(gene_tree)
# tree_nwk=gene_tree.write()

for i in gene_dic:
    if len(gene_dic[i])>1:
        # print(i)
        (gene_tree&i).up.add_child(name=gene_dic[i][1])
        (gene_tree&i).name=gene_dic[i][0]
    else:
        (gene_tree&i).name=gene_dic[i][0]
print(gene_tree.write())