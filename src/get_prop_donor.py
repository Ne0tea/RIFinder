'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-19 16:13:30
LastEditors: Ne0tea
LastEditTime: 2024-12-13 20:40:19
'''
from ete3 import Tree
import math
import src.global_variable as glv
import src.major_detect as md
import src.para_base as pb
import time

def decide_out_name(clade_subg_dic):
    if 'OUT' in clade_subg_dic:
        return 'OUT'
    elif 'Out' in clade_subg_dic:
        return 'Out'
    elif 'Outgroup' in clade_subg_dic:
        return 'Outgroup'
    elif 'outgroup' in clade_subg_dic:
        return 'outgroup'
    else:
        return None

def calculate_sister_node_clade_number_and_comapre(node:Tree,subg_clade_dic,Clade_template_dic):
    cur_clade_number={x:0 for x in Clade_template_dic}
    all_leafs=node.get_leaf_names()
    for i in all_leafs:
        if subg_clade_dic[i.split('_')[0]] not in Clade_template_dic: continue
        cur_clade_number[subg_clade_dic[i.split('_')[0]]]+=1
    ample_clade=[]
    for i in cur_clade_number:
        if cur_clade_number[i] >= Clade_template_dic[i]:
            ample_clade.append(i)
    return ample_clade
def calculate_up_node_clade_number_and_comapre(node:Tree,raw_node:Tree,subg_clade_dic,Clade_template_dic):
    cur_clade_number={x:0 for x in Clade_template_dic}
    # if not node:
    #     print(node)
    #     print(raw_node)
    raw_leafs=raw_node.get_leaf_names()
    all_leafs=node.get_leaf_names()
    all_leafs=list(set(all_leafs)-set(raw_leafs))
    for i in all_leafs:
        if subg_clade_dic[i.split('_')[0]] not in Clade_template_dic: continue
        cur_clade_number[subg_clade_dic[i.split('_')[0]]]+=1
    ample_clade=[]
    for i in cur_clade_number:
        if cur_clade_number[i] >= Clade_template_dic[i]:
            ample_clade.append(i)
    return ample_clade

def get_clade_species(node:Tree,subg_clade_dic,clade_list):
    sp_list=[]
    all_leafs=node.get_leaf_names()
    for i in all_leafs:
        if subg_clade_dic[i.split('_')[0]] in clade_list:
            sp_list.append(i.split('_')[0])
    return sp_list

def get_subfamily_list(tree:Tree):
    root = tree
    if len(root.children) == 2:
        # 获取两个分组
        group1 = root.children[0]
        group2 = root.children[1]

        # 提取两个分组的叶子节点
        subfamily1 = [leaf.name for leaf in group1.iter_leaves()]
        subfamily2 = [leaf.name for leaf in group2.iter_leaves()]
    return [subfamily1, subfamily2]

def get_subfamily_list_from_phy(temp_phylo,OUT):
    if not OUT: #if there is no OUT in sub phylogeny, just return tree
        pass
    else:
        try:
            new_root = temp_phylo&OUT
            temp_phylo.set_outgroup(new_root)
            sister_root = new_root.get_sisters()[0]
            node=temp_phylo&OUT
            node.delete()
            temp_phylo=sister_root
        except:
            return get_subfamily_list(temp_phylo)
    return get_subfamily_list(temp_phylo)

def main(gene_list:list,SG_node:Tree):
    subg_clade_dic={}
    clade_subg_dic=glv.get("clade_subg_dic")
    for key,value in clade_subg_dic.items():
        for sp in value:
            subg_clade_dic[sp]=key
    out_name=decide_out_name(clade_subg_dic)
    sub_phlo=glv.get('sub_phlo')
    subfamily_list=get_subfamily_list_from_phy(Tree(sub_phlo),out_name)
    subf1,subf2=subfamily_list

    # Clade_template_dic set the createria of clade number in subtree
    # Out name was excluded from candidate donor
    Clade_template_dic={x:int(math.sqrt(len(clade_subg_dic[x]))) for x in clade_subg_dic if x != out_name }
    # Clade_template_dic[out_name]=1

    # Every gene tree here contain outgroup gene and setted as outgroup
    gene_tree=SG_node
    gene_tree_all_clade = list(set([subg_clade_dic[x.split('_')[0]] for x in gene_tree.get_leaf_names()]))
    out_node=[x for x in gene_tree.get_leaf_names() if x.split('_')[0] in clade_subg_dic[out_name]]
    gene_tree.set_outgroup(out_node[0])
    # print(gene_tree)

    GF_asm=[x.split('_')[0] for x in gene_list]
    if len(set(GF_asm)) == 1:
        common_ancestor_sp = list(set(GF_asm))[0]
    else:
        common_ancestor_sp = '|'.join(list(set(GF_asm)))
    receive_clade = subg_clade_dic[GF_asm[0]]
    common_ancstor_node = gene_tree.get_common_ancestor(list(set(gene_list)))

    sister_node = common_ancstor_node.get_sisters()
    sister_ample_clade = []
    for i in sister_node:
        cur_sister_ample_clade = calculate_sister_node_clade_number_and_comapre(i,subg_clade_dic,Clade_template_dic)
        if len(cur_sister_ample_clade) != 0:
            sister_ample_clade = calculate_sister_node_clade_number_and_comapre(i,subg_clade_dic,Clade_template_dic)
            modi_sister_node = i

    receive_subf='subf1' if receive_clade in subf1 else 'subf2'
    sister_ample_subf=[]
    for x in sister_ample_clade :
        if x in subf1:
            sister_ample_subf.append('subf1' )
        elif x in subf2:
            sister_ample_subf.append('subf2')
    receive_sp = common_ancestor_sp
    same_subf_number = len(list(set(subfamily_list[0 if receive_clade in subf1 else 1 ]) & set(gene_tree_all_clade)))
    # print(gene_tree_all_clade, subfamily_list[0 if receive_clade in subf1 else 1 ], same_subf_number)
    cur_subf_calde=[receive_clade]
    if receive_subf in sister_ample_subf or len(sister_ample_clade) == 0:
        up_node=common_ancstor_node.up
        if not up_node:
            return None, None, None, (None,None), None
        up_node_ample_clade = calculate_up_node_clade_number_and_comapre(up_node,common_ancstor_node,subg_clade_dic,Clade_template_dic)
        up_node_ample_subf = ['subf1' if x in subf1 else 'subf2' for x in up_node_ample_clade]
        donor_clade = '|'.join(up_node_ample_clade)
        donor_sp = get_clade_species(up_node,subg_clade_dic,up_node_ample_clade)
        donor_gene = get_clade_species(up_node,subg_clade_dic,up_node_ample_clade)
        last_node=up_node
        while (len(up_node_ample_clade) == 0 or receive_subf in up_node_ample_subf) and up_node.up :

            if receive_subf in up_node_ample_subf: cur_subf_calde.extend(up_node_ample_clade)

            up_node=up_node.up
            up_node_ample_clade=calculate_up_node_clade_number_and_comapre(up_node,last_node,subg_clade_dic,Clade_template_dic)

            up_node_ample_subf=['subf1' if x in subf1 else 'subf2' for x in up_node_ample_clade]
            donor_clade='|'.join(up_node_ample_clade)
            donor_sp = get_clade_species(up_node,subg_clade_dic,up_node_ample_clade)
            donor_gene = [x for x in up_node.get_leaf_names() if subg_clade_dic[x.split('_')[0]] in up_node_ample_clade]
            last_node = up_node
    else:
        donor_clade='|'.join(sister_ample_clade)
        donor_sp = get_clade_species(modi_sister_node,subg_clade_dic,sister_ample_clade)
        donor_gene = [x for x in modi_sister_node.get_leaf_names() if subg_clade_dic[x.split('_')[0]] in sister_ample_clade]
    # print(cur_subf_calde, same_subf_number)
    if len(set(cur_subf_calde)) >= same_subf_number:
        statue='Toplogy inconsistent'
    else:
        statue='Toplogy consistent'

    if donor_sp:
        donor_sp = list(set(donor_sp))
        donor_sp = '|'.join(donor_sp)
        donor_gene = '|'.join(donor_gene)
        donor_set = (donor_sp,donor_gene)
    else:
        donor_set = (None,None)

    return statue, donor_clade, receive_clade, donor_set, receive_sp

if __name__ == "__main__":
    tree_file='E:/Bio_analysis/HGT_newcollection/4_example_recheck/Primate_test1.txt'
    # tree_file=r'E:\Bio_analysis\HGT_newcollection\1_Gene_tree_set\OG0000551.fa.rd.tml.JTT.contree'
    sup_value=70
    config = 'E:/Bio_analysis/HGT_newcollection/Primate_try/Primate_config.txt'
    # config=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    pb.main(config)
    clade_subg_dic=glv.get("clade_subg_dic")

    start_time=time.time()
    rootnested=False
    a_HGT_gene,a_HGT_gene_len,a_HGT_topnode_dic,clade_state, r_GF_dic=md.main(tree_file,sup_value,rootnested,'Pgi03A00007030',False)
    end_time=time.time()

    print(end_time - start_time)
    if clade_state:
        for i in a_HGT_gene:
            # print(AGF_top_node_dic[i])
            print(a_HGT_gene[i],len(a_HGT_gene[i]))
            # statue, donor_sp, receive_sp, donor_clade, receive_clade=main(AGF_gene_dic[i],tree_file)
            status, donor_clade, receive_clade, donor_set, receive_sp=main(a_HGT_gene[i],a_HGT_topnode_dic[i])
            print(status, donor_set[0], receive_sp, donor_clade, receive_clade)