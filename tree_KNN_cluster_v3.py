'''
Descripttion: INPUT:gene_tree NWK格式的基因树文件 
              OUTPUT:print "og_id,是否为HGT,HGTgene1 HGTgene2..."
Author: Ne0tea
version: 3.0.0
Date: 2023-06-08 22:53:43
LastEditors: Ne0tea
LastEditTime: 2023-10-21 11:55:47
'''
from sklearn.neighbors import KDTree
from sklearn.metrics import DistanceMetric
# from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import pandas as pd
from ete3 import Tree
import sys
import logging
from collections import Counter
import filter_length
import global_variable as glv

'''
根据gene_name确定科属 
INPUT:gene_name:str
OUTPUT:sub_family:str
'''

def define_sub(gene):
    clade_subg_dic=glv.get('clade_subg_dic')
    sp=gene.split("_")[0]
    for i in clade_subg_dic:
        if sp in clade_subg_dic[i]:
            return i
'''
分类Cluster中如果同时满足最大科属(dom)和最小科属(minor)比例要求,即认为发生近期HGT,并返回最小科属中的基因
INPUT:sub_list Cluster中科属
      dom 显著标准
      minor 小部分的标准
OUTPUT:statue 
       minor_gene
'''
def cluster_criteria(sub_list:list,dom:float=0.8,minor:float=0.1) -> bool:
    dic=Counter(sub_list)
    # print(dic)
    try:
        if min(dic,key=dic.get) == "OUT":
            del dic['OUT'] 
        if max(dic.values())/len(sub_list) > dom and min(dic.values())/len(sub_list) < minor:
            min_sub=min(dic.items(),key=lambda x:x[1])[0]
            return True,min_sub
        else:
            return False,None
    except ValueError:
        return False,None

'''
根据基因树中不同科属基因数量，确定聚类算法中的分类数量
'''
def cluster_define(c) -> int:
    dist = DistanceMetric.get_metric('euclidean')
    cluster_np = np.array([[6, 11, 18,25,8,2],[12, 22, 36,50,16,4]])#近单拷贝的亚科数量分布6, 11, 18,25,8,2，多拷贝OG中亚科数量分布加倍
    # print(cluster_np)
    test_np=np.append(c,cluster_np,axis=0)
    # print(test_np)
    dist_arr=dist.pairwise(test_np)[0,1:]
    if np.argmin(dist_arr) == 0:
        return 6
    else:
        return 12

'''
根据基因树,过滤长枝和串联重复基因。同时,根据leaf之间两两distance将leaf利用KNN进行分类:
INPUT: gene_tree
OUTPUT:tree_cluster_dic 树聚类结果及其中科属
       tree_cluster distance_matrix
'''
def find_HGT(gene_tree:Tree) -> dict:  
    clade_subg_dic=glv.get('clade_subg_dic')
    outbranch=filter_length.find_outlength(gene_tree)
    logging.info('{0} was delete due to extreme branch length'.format(" ".join(outbranch)))
    gene_name=[]
    sub_num={}
    data_np=np.empty((0,len(gene_tree.get_leaf_names())-len(outbranch)),float)
    for x in gene_tree:
    #针对树中长枝的过滤（V3）
        if x.name in outbranch:
            continue
        x_length=x.dist
        cur_dist=[]
        sub=define_sub(x.name)
        if sub in sub_num:
            sub_num[sub]+=1
        else:
            sub_num[sub]=1
        gene_name.append(x.name)
        for y in gene_tree:
    #针对树中长枝的过滤（V3）
            if y.name in outbranch:
                continue
            y_length=y.dist
            if x==y:
                cur_dist.append(0)
            else:
    #用于分类的distance更改为dis-node枝长
                cur_dist.append(x.get_distance(y.name)-x_length-y_length)
        # print(np.shape(cur_dist))
        cur_arry=np.array([cur_dist])
        data_np=np.append(data_np,cur_arry,axis=0)
        # gene_dist[x.name]=cur_arry

    tree_cluster = pd.DataFrame(data_np,index=gene_name)
    kdt = KDTree(data_np, leaf_size=30, metric='euclidean')
    result=kdt.query(data_np, k=3, return_distance=False)

    curSubnum=np.array([[sub_num[i] if i in sub_num else 0 for i in clade_subg_dic.keys()]])
    # print(curSubnum)

    clustering = AgglomerativeClustering(n_clusters=cluster_define(curSubnum),linkage="complete").fit(result)
    tree_cluster['cluster']=clustering.labels_
    tree_cluster['Subfamily']=[define_sub(x) for x in tree_cluster.index]
    #检查特殊分类node块
    # print(tree_cluster[['cluster','Subfamily']])
    # print(tree_cluster.loc[['Phedu_PH02Gene20425','Olat_Ola020202_2','Phedu_PH02Gene30592','Omer_OMERI01G17940_1']])

    tree_cluster_dic={}
    for i in tree_cluster['cluster'].unique():
        tree_cluster_dic["culster_"+str(i)]=[]
        for x in tree_cluster[tree_cluster['cluster'] == i].index:
            tree_cluster_dic["culster_"+str(i)].append(define_sub(x))
    return tree_cluster_dic,tree_cluster

def main(tree_file:str,dom:float=0.8,minor:float=0.1,sup_value=50):
    gene_tree=Tree(tree_file)
    go=glv.get('gene_order')
    tree_cluster_dic,tree_cluster=find_HGT(gene_tree)
    outtandem=list(filter_length.find_collinear(gene_tree,go))
    logging.info('{0} was identified as tandem genes'.format(" ".join(outtandem)))
    statue=False
    HGTgene_raw={}
    order=0
    for i in tree_cluster_dic:
        statue,min_sub=cluster_criteria(tree_cluster_dic[i],dom,minor)
        clu_num=i.split("_")[1]
        if statue:
            order+=1
            HGTgene_raw['Recent'+str(order)]=[]
            HGTgene=tree_cluster[(tree_cluster["cluster"]==int(clu_num)) & (tree_cluster["Subfamily"]==min_sub)].index
            for x in HGTgene:
                y=gene_tree&x
                while define_sub(y.get_sisters()[0].name) == define_sub(y.name):
                    y=y.up
                if y.support > sup_value:
                    HGTgene_raw['Recent'+str(order)].append(x)
    HGTgene_filt={}
    #去除HGT_list中的串联重复
    if len(HGTgene_raw) > 0:
        for i in HGTgene_raw:
            HGTgene_filt[i]=[x for x in HGTgene_raw[i] if x not in outtandem]
            if not HGTgene_filt[i]:
                del HGTgene_filt[i]
            # HGTgene_dict=[i for x in HGTgene_raw if i not in outtandem]
        return HGTgene_filt,True
        # print(og,True,gene_line)
    else:
        return None,False
        # print(og,False,"NO HGT")


if __name__ == '__main__':
    tree_file=sys.argv[1]
    dom=sys.argv[2]
    minor=sys.argv[3]
    sup_value=sys.argv[4]
    tree_file=r"D:\study resource\HGT\4_significant_gene\OG0003936_ali.mafft_rmdup_trimal_JTT.treefile"
    main(tree_file,dom,minor,sup_value)


# pd.options.display.width=10000
# pd.options.display.max_rows = None
# # pd.options.display.max_colwidth=10
# pd.options.display.precision=3
# pd.options.display.max_columns = None