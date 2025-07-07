'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-16 17:11:48
LastEditors: Ne0tea
LastEditTime: 2024-09-16 17:12:43
'''
import pandas as pd
from ete3 import Tree
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