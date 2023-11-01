'''
cripttion: 
Author: Ne0tea
version: 
Date: 2023-05-16 21:53:43
LastEditors: Ne0tea
LastEditTime: 2023-10-21 11:42:28
'''
import sys
from ete3 import Tree
import re
import copy
import global_variable as glv

def define_sub(gene):
    clade_subg_dic=glv.get('clade_subg_dic')
    sp=gene.split("_")[0]
    for i in clade_subg_dic:
        if sp in clade_subg_dic[i]:
            return i

##gene tree collapse criteria 
def collapsed_leaf(node):
    sub_list=[]
    for i in node2labels1[node]:
        # print(i.name)
        if i != "":
            sub=define_sub(i)
            sub_list.append(sub)
        else:
            sub_list.append("")

    sub_list=list(filter(None,sub_list))
    # print(len(set(sub_list)))
    if len(set(sub_list)) == 1:
        node.name=sub_list[0]
        return True
    else:
        return False

def collapsed_leaf_round2(node):
    if len(node2labels2[node]) == 1:
        node.name=list(node2labels2[node])[0]
        return True
    else:
        return False




count=0
def replace(m):   # 在函数中使用外部定义的计数器
    global count
    count += 1     # 计数器加1
    return ')' + 'node' + str(count) + ':'   # 替换为新节点名称

def main(tree_file:str,standard:int=50):
    gene_tree=Tree(tree_file)
    gene_tree2=open(tree_file,'r')
    gene_tree2=gene_tree2.readline()
    # print(gene_tree2)
    pattern = re.compile(r'\)(\d*):')
    inter_nametree=re.sub(pattern,replace, gene_tree2)
    new_tree=Tree(inter_nametree,format=1)
    # print(new_tree)
    global node2labels1
    node2labels1 = new_tree.get_cached_content(store_attr="name")
    collapse_tree = Tree(new_tree.write(is_leaf_fn=collapsed_leaf))
    round1_leaf_num=len(collapse_tree.get_leaves())
    # collapse_tree.write(format=1, outfile=tree_file+"_collapse")
    try:
        out=collapse_tree&"OUT"
        collapse_tree.set_outgroup(out)
    except:
        pass
    # print(collapse_tree.write())
    collapse_tree_intername=re.sub(pattern,replace, collapse_tree.write())
    new_tree=Tree(collapse_tree_intername,format=1)
    global node2labels2
    node2labels2 = new_tree.get_cached_content(store_attr="name")
    # print(node2labels)
    collapse_tree2 = Tree(new_tree.write(is_leaf_fn=collapsed_leaf_round2))
    topo_num=0
    new_tree=copy.copy(collapse_tree2)
    for topo in collapse_tree2.traverse("postorder"):
        topo_num+=1
        if len(topo.get_leaf_names())==2 and topo.children[0].name:
            layer=0
            topo2=topo
            while topo2:
                if layer==0:
                    sister_node=topo2.get_sisters()
                    # print(sister_node,"|",topo2.get_leaf_names())
                    if sister_node and sister_node[0].name in topo2.get_leaf_names():
                        # print(sister_node[0].name)
                        topo2.remove_sister(sister_node[0])
                        # print(collapse_tree2)
                    last=sister_node[0].name
                    topo2=topo2.up
                    # topo2=topo2.up
                    layer+=1
                    continue
                if topo2.get_sisters():
                    sister_node=topo2.get_sisters()
                    if sister_node[0].name == last and sister_node[0].name:
                        topo2.remove_sister(sister_node[0])
                    last=sister_node[0].name
                    topo2=topo2.up 
                    layer+=1
                else:
                    topo2=topo2.up
    round2_leaf_num=len(collapse_tree2.get_leaves())
    # print(tree_file,topo_num,round1_leaf_num,round2_leaf_num)
    if round2_leaf_num > standard:
        return round2_leaf_num,False
    else:
        return round2_leaf_num,True

if __name__ == "__main__":
    # tree_file=sys.argv[1]
    # tf_file=sys.argv[2]
    tree_file=r"D:\study resource\HGT\4_significant_gene\OG0003936_ali.mafft_rmdup_trimal_JTT.treefile"
    # tf_file=r"D:\Bio_analysis\Grass_horizatal_transversion\test\costume_check\OG0014619_transfer_pair"
    main(tree_file)