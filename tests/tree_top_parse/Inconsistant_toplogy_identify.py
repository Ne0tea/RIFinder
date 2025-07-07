'''
Descripttion: Determining the events of OG-like structures
Author: Ne0tea
version: 
Date: 2024-09-12 10:01:15
LastEditors: Ne0tea
LastEditTime: 2024-11-02 11:41:34
'''

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

'''
contary_step: indicating hierarchical level
TGF: Trickle gene flow
AGF: Ancient gene flow
'''
def rooted_toplogy(process,subfam_list,out_name:str,cur_clade_process):
    subf1,subf2=subfam_list
    contary_switch=False
    contary_count=0
    contary_step=0
    last_value=[]
    plus_not_None=0
    for index,value in enumerate(process):
        plus_value=value[len(last_value):]
        if not plus_value: # Skip empty plus
            continue
        plus_not_None+=1
        if plus_not_None==2:
            return_index=index-1
            for back_clade in cur_clade_process[:index-1][::-1]:
                if back_clade != cur_clade_process[0]:
                    return_index-=1
                else:
                    break

        cur_process_clade=['subf1' if x in subf1 else 'subf2' for x in plus_value if x != out_name]
        last_value=value
        if plus_not_None == 1:
            cur_clade='subf1' if plus_value[0] in subf1 else 'subf2'
            contary_clade='subf2' if plus_value[0] in subf1 else 'subf1'
            contary_clade_count = len(subf1) if contary_clade == 'subf1' else len(subf2)

            #Use Cayley for calculate candiate top number
            contary_clade_top_count=contary_clade_count**(contary_clade_count-2)
            continue
        if contary_clade in cur_process_clade and not contary_switch:
            if plus_not_None != 2:
                contary_switch = True
                # print(contary_clade_top_count, cur_process_clade.count(contary_clade) / len(cur_process_clade))
                contary_step += contary_clade_top_count / cur_process_clade.count(contary_clade) / len(cur_process_clade)
                if cur_clade in cur_process_clade:
                    break
            else:
                #less value for sister node in the sg tree
                if cur_clade not in cur_process_clade:
                    contary_step += 1 / 2 / cur_process_clade.count(contary_clade) / len(cur_process_clade)
            contary_count+=cur_process_clade.count(contary_clade)
            continue
        # print(contary_step)
        if contary_switch:
            if contary_clade in cur_process_clade:
                contary_step += contary_clade_top_count / cur_process_clade.count(contary_clade) / len(cur_process_clade)
            if cur_clade in cur_process_clade:
                break
            else:
                contary_count+=cur_process_clade.count(contary_clade)

    # print(contary_count, contary_step, contary_clade_top_count)
    if contary_count >= contary_clade_count:
        if contary_step < contary_clade_top_count:
            return 'TGF', cur_clade, return_index
        else:
            return 'AGF', contary_clade, return_index
    else:
        return 'Inconstant_clade', None, None

def decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name,cur_clade_process):
    first_process= [xx for xx in cur_ample_process if xx]
    if not first_process:
        return 'Outlier_clade', None, None
    first_process = first_process[0]
    if len(first_process) != 1:
        return 'Outlier_clade', None, None
    # print(cur_ample_process)
    if out_name in cur_ample_process[-1]:
        toplogy,clade_name,clade_index = rooted_toplogy(cur_ample_process,subfamily_list,out_name,cur_clade_process)
        # print(clade_index)
    else:
        # toplogy,clade_name,clade_index=unroot_toplogy(cur_ample_process,subfamily_list,out_name)
        return 'Norooted', None, None
    return toplogy,clade_name,clade_index


if __name__ == "__main__":
    tree_file=r'E:\Bio_analysis\HGT_newcollection\1_Gene_tree_set\OG0000551.fa.rd.tml.JTT.contree'
    sup_value=70
    gene_order=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\All_gene.sord'
    all_process=[]
    if 1:
        # INconstant_clade
        process1=[['ORY'], ['ORY', 'BAM'], ['ORY', 'BAM', 'CHL'], ['ORY', 'BAM', 'CHL', 'POO', 'PAN'], ['PAN', 'CHL', 'POO', 'ORY', 'BAM', 'OUT']]
        # INconstant_clade
        process2=[['CHL'],['CHL','ORY','BAM'], ['CHL','ORY','BAM', 'POO', 'PAN'],['CHL','ORY','BAM', 'POO', 'PAN', 'OUT']]
        # INconstant_clade
        process3=[['PAN'],['PAN','POO'],['PAN','POO','PAN'],['PAN','POO','PAN','CHL','ORY','BAM'],['PAN','POO','PAN','CHL','ORY','BAM'],['PAN','POO','PAN','CHL','ORY','BAM','OUT']]
        # AGF
        process4=[['POO'],['POO','PAN'],['POO','PAN','PAN'],['POO','PAN','PAN','CHL','ORY','BAM'],['POO','PAN','PAN','CHL','ORY','BAM'],['POO','PAN','PAN','CHL','ORY','BAM','OUT']]
        # TGF
        process5=[['ORY'], ['ORY', 'BAM'], ['ORY', 'BAM', 'POO'], ['ORY', 'BAM', 'POO', 'CHL', 'PAN'], ['ORY', 'BAM', 'POO', 'CHL', 'PAN', 'OUT']]
        # TGF
        process6=[['CHL'],['CHL','PAN'], ['CHL','PAN','ORY','BAM', 'POO'],['CHL','PAN','ORY','BAM', 'POO', 'OUT']]
        # AGF
        process7=[['PAN'],['PAN','CHL'],['PAN','CHL','ORY'],['PAN','CHL','ORY','BAM'],['PAN','CHL','ORY','BAM','POO'],['PAN','CHL','ORY','BAM','POO','OUT']]
        # AGF
        process8=[['PAN'],['PAN','CHL'],['PAN','CHL','ORY'],['PAN','CHL','ORY','BAM','POO'],['PAN','CHL','ORY','BAM','POO'],['PAN','CHL','ORY','BAM','POO','OUT']]
        process9=[['PAN'],['PAN','CHL'],['PAN','CHL','ORY'],['PAN','CHL','ORY','POO'],['PAN','CHL','ORY','POO','BAM'],['PAN','CHL','ORY','POO','BAM','OUT']]
        process10=[['ORY'],['ORY','POO'],['ORY','POO','BAM'],['ORY','POO','BAM','PAN'],['ORY','POO','BAM','PAN','CHL'],['ORY','POO','BAM','PAN','CHL','OUT']]
        # AGF
        process11=[['POO'],['POO','PAN'],['POO','PAN','CHL'],['POO','PAN','CHL','ORY','BAM'],['POO','PAN','CHL','ORY','BAM','OUT']]
        # INconstant_clade
        process12=[['ORY'],['ORY','BAM'],['ORY','BAM','PAN'],['ORY','BAM','PAN','POO','CHL'],['ORY','BAM','PAN','POO','CHL','OUT']]
        # Inconstant_clade
        process13=[['BAM'],['BAM','ORY'],['BAM','ORY','PAN'],['BAM','ORY','PAN'],['BAM','ORY','PAN','PAN','POO','CHL'],['BAM','ORY','PAN','PAN','POO','CHL','OUT']]
        # Inconstant_clade
        process14=[['PAN'], ['PAN', 'PAN'], ['PAN', 'PAN', 'CHL', 'POO', 'PAN'], ['PAN', 'PAN', 'CHL', 'POO', 'PAN', 'POO', 'BAM', 'ORY'], ['PAN', 'PAN', 'CHL', 'POO', 'PAN', 'POO', 'BAM', 'ORY', 'OUT']]
        # AGF
        process15=[['POO'], ['POO'], ['POO', 'PAN'], ['POO', 'PAN'], ['POO', 'PAN', 'PAN'], ['POO', 'PAN', 'PAN'], ['POO', 'PAN', 'PAN', 'BAM', 'ORY', 'CHL'], ['POO', 'PAN', 'PAN', 'BAM', 'ORY', 'CHL', 'OUT']]
        process16=[['PAN'], ['PAN','CHL'], ['PAN','CHL','BAM','ORY','POO'], ['PAN','CHL','BAM','ORY','POO','BAM'],  ['PAN','CHL','BAM','ORY','POO','BAM','OUT']]
        # AGF
        process17=[['POO'],['POO','PAN'],['POO','PAN','PAN'],['POO','PAN','PAN','CHL'],['POO','PAN','PAN','CHL','ORY','BAM'],['POO','PAN','PAN','CHL','ORY','BAM','OUT']]

    if 1:
        cur_clade_process1 = [['ORY'], ['BAM'], ['CHL'], ['POO', 'PAN'], ['OUT']]
        cur_clade_process2 = [['CHL'], ['ORY', 'BAM'], ['POO', 'PAN'], ['OUT']]
        cur_clade_process3 = [['PAN'], ['POO'], ['PAN'], ['CHL', 'ORY', 'BAM'], [], ['OUT']]
        cur_clade_process4 = [['POO'], ['PAN'], ['PAN'], ['CHL', 'ORY', 'BAM'], [], ['OUT']]
        cur_clade_process5 = [['ORY'], ['BAM'], ['POO'], ['CHL', 'PAN'], ['OUT']]
        cur_clade_process6 = [['CHL'], ['PAN'], ['ORY', 'BAM', 'POO'], ['OUT']]
        cur_clade_process7 = [['PAN'], ['CHL'], ['ORY'], ['BAM'], ['POO'], ['OUT']]
        cur_clade_process8 = [['PAN'], ['CHL'], ['ORY'], ['BAM', 'POO'], [], ['OUT']]
        cur_clade_process9 = [['PAN'], ['CHL'], ['ORY'], ['POO'], ['BAM'], ['OUT']]
        cur_clade_process10 = [['ORY'], ['POO'], ['BAM'], ['PAN'], ['CHL'], ['OUT']]
        cur_clade_process11 = [['POO'], ['PAN'], ['CHL'], ['ORY', 'BAM'], ['OUT']]
        cur_clade_process12 = [['ORY'], ['BAM'], ['PAN'], ['POO', 'CHL'], ['OUT']]
        cur_clade_process13 = [['BAM'], ['ORY'], ['PAN'], [], ['PAN', 'POO', 'CHL'], ['OUT']]
        cur_clade_process14 = [['PAN'], ['PAN'], ['CHL', 'POO', 'PAN'], ['POO', 'BAM', 'ORY'], ['OUT']]
        cur_clade_process15 = [['POO'], [], ['PAN'], [], ['PAN'], [], ['BAM', 'ORY', 'CHL'], ['OUT']]
        cur_clade_process16 = [['PAN'], ['CHL'], ['BAM', 'ORY', 'POO'], ['BAM'], ['OUT']]
        cur_clade_process17 = [['POO'], ['PAN'], ['PAN'], ['CHL'], ['ORY', 'BAM'], ['OUT']]
    cur_clade_process=[['POO'], ['PAN'], ['POO'], ['POO'], ['POO'], ['POO'], ['BAM'], ['PAN', 'CHL'], ['ORY'], ['ORY'], ['ORY'], ['ORY'], ['POO', 'ORY'], ['OUT']]
    # process=[['PAN'], ['PAN', 'CHL'], ['PAN', 'CHL', 'PAN', 'CHL'], ['PAN', 'CHL', 'PAN', 'CHL', 'POO', 'BAM', 'ORY'], ['PAN', 'CHL', 'PAN', 'CHL', 'POO', 'BAM', 'ORY', 'POO'], ['PAN', 'CHL', 'PAN', 'CHL', 'POO', 'BAM', 'ORY', 'POO', 'POO', 'PAN', 'CHL', 'BAM'], ['PAN', 'CHL', 'PAN', 'CHL', 'POO', 'BAM', 'ORY', 'POO', 'POO', 'PAN', 'CHL', 'BAM', 'ORY'], ['PAN', 'CHL', 'PAN', 'CHL', 'POO', 'BAM', 'ORY', 'POO', 'POO', 'PAN', 'CHL', 'BAM', 'ORY', 'OUT']]
    process=[['POO'], ['POO', 'PAN'], ['POO', 'PAN', 'POO'], ['POO', 'PAN', 'POO', 'POO'], ['POO', 'PAN', 'POO', 'POO'], ['POO', 'PAN', 'POO', 'POO', 'POO'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL', 'ORY'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL', 'ORY'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL', 'ORY', 'POO'], ['POO', 'PAN', 'POO', 'POO', 'POO', 'BAM', 'PAN', 'CHL', 'ORY', 'POO', 'OUT']]
    all_process=[process1,process2,process3,process4,process5,process6,process7,process8,process9,process10,process11,process12,process13,process14,process15,process16,process17]
    all_clade_process=[cur_clade_process1,cur_clade_process2,cur_clade_process3,cur_clade_process4,cur_clade_process5,cur_clade_process6,cur_clade_process7,cur_clade_process8,cur_clade_process9,cur_clade_process10,cur_clade_process11,cur_clade_process12,cur_clade_process13,cur_clade_process14,cur_clade_process15,cur_clade_process16,cur_clade_process17]
    # all_process=[process]
    subfamily_list=[['CHL', 'PAN'],['POO', 'ORY', 'BAM']]
    clade_template_dic={'BAM': 5, 'ORY': 3, 'POO': 5, 'PAN': 6, 'CHL': 3, 'OUT':1}
    count=0
    for index,process in enumerate(all_process):
        count+=1
        xx=decide_toplogy_from_process(process,subfamily_list,'OUT',all_clade_process[index])
        print(xx)
        print('process',count)