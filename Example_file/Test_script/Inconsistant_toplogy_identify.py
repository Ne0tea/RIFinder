'''
Descripttion: Determining the events of OG-like structures
Author: Ne0tea
version: 
Date: 2024-09-12 10:01:15
LastEditors: Ne0tea
LastEditTime: 2024-09-14 16:45:01
'''
from itertools import combinations
import para_base  as pb
import global_variable as glv
from ete3 import Tree
import numpy as np

def matrix_in_list(matrix, matrix_list):
    """判断矩阵是否在矩阵组成的列表中"""
    return any(np.array_equal(matrix, mat) for mat in matrix_list)

def matrix_determine(process_matrix):
    for i, row in enumerate(process_matrix):
        if i == 0 and row != [0,1] and row != [1,0]:
            raise 'It seems like a problem'
        if row[0] > 0 and row[1] > 0:
            return 'Outlier_clade',None
        pass

def decide_out_name(clade_subg_dic):
    if 'OUT' in clade_subg_dic:
        return 'OUT'
    elif 'Out' in clade_subg_dic:
        return 'Out'
    elif 'Outgroup' in clade_subg_dic:
        return 'Outgroup'
    elif 'outgroup' in clade_subg_dic:
        return 'outgroup'

def fragmented_toplogy(process,subfam_list):
    for step in process:
        print(step)
    pass

def unroot_toplogy(process,subfam_list,out_name):
    subf1,subf2=subfam_list
    contary_count=0
    contray_step=0
    last_value=[]
    plus_not_None=0
    for index,value in enumerate(process):
        plus_value=list(set(value) - set(last_value))
        if not plus_value: # Skip empty plus
            continue
        plus_not_None+=1
        if plus_not_None==2:
            return_index=index-1
        cur_process_clade = ['subf1' if x in subf1 else 'subf2' for x in plus_value if x != out_name]
        last_value = value
        if index == 0:
            cur_clade = 'subf1' if plus_value[0] in subf1 else 'subf2'
            cur_clade_count=len(subf1) if cur_clade == 'subf1' else len(subf2)
            contray_clade = 'subf2' if plus_value[0] in subf1 else 'subf1'
            contary_clade_count=len(subf1) if contray_clade == 'subf1' else len(subf2)
            same_count=1
            same_step=1
            continue
        # print(cur_process_clade)
        if len(set(cur_process_clade)) == 2:
            break
        else:
            if cur_clade in cur_process_clade:
                same_count += cur_process_clade.count(cur_clade)
                same_step += 1
            else:
                contary_count += cur_process_clade.count(contray_clade)
                contray_step+=1
                same_forze=same_step
    if contary_count == contary_clade_count:
        if contary_clade_count >= 3 and contray_step >= 2 and same_forze == 1:
            return 'AGF', contray_clade, return_index
        else:
            return 'TGF', cur_clade, return_index
    elif same_count == cur_clade_count and contary_count < contary_clade_count:
         return 'TGF', cur_clade, return_index
    else:
        return 'Inconstant_clade', None, None


'''
contray_step: indicating hierarchical level
TGF: Trickle gene flow
AGF: Ancient gene flow
'''
def rooted_toplogy(process,subfam_list,out_name):
    subf1,subf2=subfam_list
    contary_switch=False
    contary_count=0
    contray_step=0
    last_value=[]
    plus_not_None=0
    for index,value in enumerate(process):
        plus_value=list(set(value) - set(last_value))
        if not plus_value: # Skip empty plus
            continue
        plus_not_None+=1
        if plus_not_None==2:
            return_index=index
        cur_process_clade=['subf1' if x in subf1 else 'subf2' for x in plus_value if x != out_name]
        last_value=value
        if index == 0:
            cur_clade='subf1' if plus_value[0] in subf1 else 'subf2'
            contray_clade='subf2' if plus_value[0] in subf1 else 'subf1'
            continue
        if contray_clade in cur_process_clade and not contary_switch:
            contary_switch=True
            contray_step+=1
            contary_count+=cur_process_clade.count(contray_clade)
            continue
        if contary_switch:
            contray_step+=1
            if cur_clade in cur_process_clade:
                contary_count+=cur_process_clade.count(contray_clade)
                # contary_index=index-1
                break
            else:
                contary_count+=cur_process_clade.count(contray_clade)
                # contary_index=index
    contary_clade_count=len(subf1) if contray_clade == 'subf1' else len(subf2)
    # print(contary_count,contray_step)
    if contary_count == contary_clade_count:
        if contray_step <= 2:
            return 'TGF', cur_clade, return_index
        else:
            return 'AGF', contray_clade, return_index
    else:
        return 'Inconstant_clade', None, None

def decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name,Clade_template_dic):
    first_process= [xx for xx in cur_ample_process if xx][0]
    if len(first_process) != 1:
        return 'Outlier_clade', None, None

    # OUT nested in subphy wasn't under consideration
    if out_name in cur_ample_process[-1] and out_name not in cur_ample_process[-2]:
        toplogy,clade_name,clade_index=rooted_toplogy(cur_ample_process,subfamily_list,out_name)
    elif out_name in cur_ample_process[-2]:
        toplogy,clade_name,clade_index=rooted_toplogy(cur_ample_process,subfamily_list,out_name)
    # elif any(item not in cur_ample_process[-1] for item in Clade_template_dic):
    #     toplogy,clade_name,clade_index=fragmented_toplogy(cur_ample_process,subfamily_list,out_name)
    else:
        toplogy,clade_name,clade_index=unroot_toplogy(cur_ample_process,subfamily_list,out_name)
    return toplogy,clade_name,clade_index


if __name__ == "__main__":
    tree_file=r'E:\Bio_analysis\HGT_newcollection\1_Gene_tree_set\OG0000551.fa.rd.tml.JTT.contree'
    sup_value=50
    config=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    gene_order=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\All_gene.ord'
    pb.main(config,gene_order)
    process=[['CHL'], ['CHL', 'ORY'], ['CHL', 'BAM', 'ORY','POO'], ['CHL', 'PAN', 'ORY', 'BAM', 'POO']]
    subfamily_list=[['CHL', 'PAN'],['POO', 'ORY', 'BAM']]
    clade_template_dic={'BAM': 7, 'ORY': 3, 'POO': 7, 'PAN': 8, 'CHL': 2, 'OUT': 1}
    xx=decide_toplogy_from_process(process,subfamily_list,'OUT',clade_template_dic)
    print(xx)