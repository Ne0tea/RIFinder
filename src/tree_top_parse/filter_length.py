'''
Descripttion: 对传递进来的物种树(ete3)对象，识别长枝及串联重复。
Author: Ne0tea
version: 
Date: 2023-06-25 20:27:53
LastEditors: Ne0tea
LastEditTime: 2024-10-13 23:06:18
'''
from ete3 import Tree
from scipy import stats
import numpy as np
import re

def distribution_check(region:tuple,l:int) -> bool:
    # print(region)
    if l > region[1] or l < region[0]:
        return True
def find_outlength(gene_tree:Tree,outsp=[]) -> list:
    name_length={}
    outliner=[]
    length_set1=np.empty(shape=(0,1))
    length_set2=np.empty(shape=(0,1))
    length_set3=np.empty(shape=(0,1))
    for i in gene_tree:
        upper1=i.up
        try:
            upper2=upper1.up
            upper3=upper2.up
            upper_check=upper3.up
        except AttributeError:
            outliner.append(i.name)
            continue
        dis1=upper1.get_distance(i.name)
        dis2=upper2.get_distance(i.name)
        dis3=upper3.get_distance(i.name)
        # dis=gene_tree.get_distance(i.name)
        length_set1=np.append(length_set1, values=dis1)
        length_set2=np.append(length_set2, values=dis2)
        length_set3=np.append(length_set3, values=dis3)
        name_length[i.name]=[dis1,dis2,dis3]

    mean1, std1 = length_set1.mean(), length_set1.std(ddof=1)
    mean2, std2 = length_set2.mean(), length_set2.std(ddof=1)
    mean3, std3 = length_set3.mean(), length_set3.std(ddof=1)
    conf_intveral1 = stats.norm.interval(0.99, loc=mean1, scale=std1)
    conf_intveral2 = stats.norm.interval(0.99, loc=mean2, scale=std2)
    conf_intveral3 = stats.norm.interval(0.99, loc=mean3, scale=std3)
    for x in name_length:
        if distribution_check(conf_intveral1,name_length[x][0]) or \
            distribution_check(conf_intveral2,name_length[x][1]) or \
            distribution_check(conf_intveral3,name_length[x][2]):
            outliner.append(x)
    outliner = [x for x in outliner if x.split('_')[0] not in outsp]
    return outliner
def find_longest_numeric_sequence(string):
    pattern = r'(?<![Cc]hr)\d+'  # 匹配一个或多个数字
    matches = re.findall(pattern, string)
    longest_sequence = max(matches, key=len) if matches else ''
    return longest_sequence

if __name__ == '__main__':
    tree_file=r"D:\study resource\HGT\4_significant_gene\OG0003195_ali.mafft_rmdup_trimal_JTT.treefile"
    order_file=r"E:\Bio_analysis\Grass_horizatal_transversion\subg_gene.order"
    gene_tree=Tree(tree_file)
    # outbranch=find_outlength(gene_tree)
    # outtandem=find_collinear(gene_tree,order_file)
    # print(outtandem)