'''
Descripttion: 对传递进来的物种树(ete3)对象，识别长枝及串联重复。
Author: Ne0tea
version: 
Date: 2023-06-25 20:27:53
LastEditors: Ne0tea
LastEditTime: 2023-10-20 17:11:52
'''
import pandas as pd
from ete3 import Tree
from scipy import stats
import numpy as np
import re

def distribution_check(region:tuple,l:int) -> bool:
    # print(region)
    if l > region[1] or l < region[0]:
        return True
def find_outlength(gene_tree:Tree) -> list:
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
        # if i.name=="Atau_AET5Gv21056700" or i.name=="Pten_Chr0600824":
        #     print(name_length[i.name])
    mean1, std1 = length_set1.mean(), length_set1.std(ddof=1)
    mean2, std2 = length_set2.mean(), length_set2.std(ddof=1)
    mean3, std3 = length_set3.mean(), length_set3.std(ddof=1)
    # print(length_set)
    conf_intveral1 = stats.norm.interval(0.99, loc=mean1, scale=std1)
    # print(conf_intveral1)
    conf_intveral2 = stats.norm.interval(0.99, loc=mean2, scale=std2)
    # print(conf_intveral2)
    conf_intveral3 = stats.norm.interval(0.99, loc=mean3, scale=std3)
    # print(conf_intveral3)
    for x in name_length:
        if distribution_check(conf_intveral1,name_length[x][0]) or \
            distribution_check(conf_intveral2,name_length[x][1]) or \
            distribution_check(conf_intveral3,name_length[x][2]):
            outliner.append(x)
    return outliner
def find_longest_numeric_sequence(string):
    pattern = r'(?<![Cc]hr)\d+'  # 匹配一个或多个数字
    matches = re.findall(pattern, string)
    longest_sequence = max(matches, key=len) if matches else ''
    return longest_sequence
def find_collinear(gene_tree:Tree,gene_order:pd) -> list:
    order_data=pd.read_table(gene_order,names=['gene','Chr','order'])
    # tandem_set=[]
    leaf_name=gene_tree.get_leaf_names()
    leaf_loc=order_data.loc[order_data['gene'].isin(leaf_name)]
    # print(leaf_loc)
    tandem_target=[]
    for i in set(leaf_loc['Chr']):
        c_loc=leaf_loc.loc[leaf_loc['Chr']==i,'order'].tolist()
        c_loc.sort()
        # print(c_loc)
        if len(c_loc)>=3:
            result = list(filter(lambda x: c_loc[x+1] - c_loc[x] < 10 or c_loc[x] - c_loc[x-1] < 10, range(1, len(c_loc)-1)))
            if c_loc[1] - c_loc[0] < 10:
                result.insert(0, 0)  # 添加第一个元素
            if c_loc[-1] - c_loc[-2] < 10:
                result.append(len(c_loc)-1)  # 添加最后一个元素
        elif len(c_loc)==2:
            if c_loc[1] - c_loc[0] < 10:
                result=[0,1]
            else:
                continue
        else:
            continue
        c_result=[c_loc[x] for x in result]
        # print(leaf_loc.loc[(leaf_loc['Chr']==i)&(leaf_loc['order'].isin(c_result))]['gene'].tolist())
        if tandem_target:
            tandem_target.extend(leaf_loc.loc[(leaf_loc['Chr']==i)&(leaf_loc['order'].isin(c_result))]['gene'].tolist())
        else:
            tandem_target=leaf_loc.loc[(leaf_loc['Chr']==i)&(leaf_loc['order'].isin(c_result))]['gene'].tolist()
    return list(set(tandem_target))


if __name__ == '__main__':
    tree_file=r"D:\study resource\HGT\4_significant_gene\OG0003195_ali.mafft_rmdup_trimal_JTT.treefile"
    order_file=r"E:\Bio_analysis\Grass_horizatal_transversion\subg_gene.order"
    gene_tree=Tree(tree_file)
    # outbranch=find_outlength(gene_tree)
    outtandem=find_collinear(gene_tree,order_file)
    print(outtandem)