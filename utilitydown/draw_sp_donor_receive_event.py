'''
Descripttion: generate gene flow gene cluster and named as AGFD_stat
Author: Ne0tea
version: 
Date: 2024-09-29 21:50:55
LastEditors: Ne0tea
LastEditTime: 2024-11-13 00:36:26
'''
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator

# template_color={'PLA':'#e9c46a','LOR':'#264653','HOM':'#e76f51','CHL':'#2a9d8f','CER':'#f4a261'}
template_color={"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1"}
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
out_file=r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out.consistance.significant.all.out'
config_file = r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
clade_subg_dic={}
subg_clade_dic={}
with open(config_file,'r') as c_file:
    phlo=0
    subg=0
    for i in c_file:
        if i.startswith('>'):
            if i.startswith('>phylogeny'):
                phlo,subg=1,0
            if i.startswith('>subg list'):
                phlo,subg=0,1
        else:
            if phlo:
                sub_phlo=i.strip()
            if subg:
                line=i.strip().split()
                clade_subg_dic[line[0]]=line[1:]
                for subg in line[1:]:
                    subg_clade_dic[subg]=line[0]

sp_ord=[]
for x in clade_subg_dic:
    if x != 'OUT':
        sp_ord.extend(clade_subg_dic[x])

# if 'calculate based on event':
#     sp_receive_dic, sp_donor_dic = {}, {}
#     with open(out_file, 'r') as of:
#         for i in of:
#             cur_line=i.strip().split('\t')
#             '''
#             out file format example:
#             donor_clade, receive_clade, donor_sp, receive_sp
#             GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
#             '''
#             gf_id, donor_clade, receive_clade, donor_sp, donor_gene, receive_sp, receive_gene = cur_line[0], cur_line[4],cur_line[5],cur_line[6],cur_line[7],cur_line[8],cur_line[9]
#             for sp in donor_sp.split('|'):
#                 if sp in sp_receive_dic:
#                     for clade in receive_clade.split('|'):
#                         if clade in sp_receive_dic[sp]:
#                             sp_receive_dic[sp][clade] += 1
#                         else:
#                             sp_receive_dic[sp][clade] = 1
#                 else:
#                     for clade in receive_clade.split('|'):
#                         if sp in sp_receive_dic:
#                             sp_receive_dic[sp][clade] = 1
#                         else:
#                             sp_receive_dic[sp] = {clade:1}
#             for sp in receive_sp.split('|'):
#                 if sp in sp_donor_dic:
#                     for clade in donor_clade.split('|'):
#                         if clade in sp_donor_dic[sp]:
#                             sp_donor_dic[sp][clade] += 1
#                         else:
#                             sp_donor_dic[sp][clade] = 1
#                 else:
#                     for clade in donor_clade.split('|'):
#                         if sp in sp_donor_dic:
#                             sp_donor_dic[sp][clade] = 1
#                         else:
#                             sp_donor_dic[sp] = {clade:1}

#     sp_receive_df = pd.DataFrame.from_dict(sp_donor_dic, orient='index')
#     sp_donor_df = pd.DataFrame.from_dict(sp_receive_dic, orient='index')
#     sp_receive_df['Clade'] = sp_receive_df.index.map(subg_clade_dic)
#     sp_donor_df['Clade'] = sp_donor_df.index.map(subg_clade_dic)
#     sp_receive_df = sp_receive_df.reindex(sp_ord)
#     sp_donor_df = sp_donor_df.reindex(sp_ord)
#     sp_receive_df = sp_receive_df.fillna(0)
#     sp_donor_df = sp_donor_df.fillna(0)

#     plt.figure(figsize=(6, 16))
#     sns.set_style('whitegrid')

#     group_gap = 1.5
#     y_positions = []
#     current_position = 0
#     for group in sp_receive_df['Clade'].unique():
#         group_data = sp_receive_df[sp_receive_df['Clade'] == group]
#         y_positions.extend(np.arange(current_position, current_position + len(group_data)))
#         current_position += len(group_data) + group_gap

#     last_value = [0] * len(sp_receive_df)
#     for group in sp_receive_df['Clade'].unique():
#         plt.barh(y_positions, sp_receive_df[group].to_list(), left=last_value, linewidth=0, 
#                 color = template_color[group], edgecolor='none')
#         last_value = [x+last_value[idx] for idx,x in enumerate(sp_receive_df[group].to_list())]

#     last_value = [0] * len(sp_receive_df)
#     for group in sp_donor_df['Clade'].unique():
#         plt.barh(y_positions, [-num for num in sp_donor_df[group].to_list()], left=last_value, linewidth=0, 
#                 color=template_color[group], label=group, edgecolor='none')
#         last_value = [-num+last_value[idx] for idx,num in enumerate(sp_donor_df[group].to_list())]


#     xmin_value = sp_donor_df[list(template_color.keys())].sum(axis=1).max()
#     xmax_value = sp_receive_df[list(template_color.keys())].sum(axis=1).max()
#     ax = plt.gca()
#     current_position = -0.5
#     anno_rectangle = 20
#     for group in sp_receive_df['Clade'].unique():
#         group_data = sp_receive_df[sp_receive_df['Clade'] == group]
#         width = len(group_data)
#         rect = patches.Rectangle((-xmin_value-anno_rectangle, current_position), anno_rectangle, width,
#                                 linewidth=0.5, edgecolor='black', facecolor=template_color[group], alpha=0.6)
#         ax.add_patch(rect)
#         ax.text(-xmin_value-anno_rectangle/2, current_position+width/2, group,
#                 rotation=-90, horizontalalignment='center', verticalalignment='center', fontweight='bold')
#         current_position += len(group_data) + group_gap

#     plt.yticks(y_positions, list(sp_receive_df.index.values))
#     # plt.ylabel('Value', labelpad=15)
#     ax.set_xlim(-xmin_value-20, xmax_value)  # 限制 x 轴范围
#     ax.set_ylim(-1, 128)
#     plt.axvline(0, color='white', linewidth=1)
#     plt.legend()
#     plt.gca().xaxis.set_major_locator(MaxNLocator(10))
#     sns.despine(left=True, bottom=True)
#     plt.legend(loc='upper right')
#     plt.tight_layout()
#     # plt.show()
#     plt.savefig(r'E:\Bio_analysis\HGT_newcollection\3_outfile_stat\AGFD_event_stat.pdf', bbox_inches='tight', format='pdf')
#     plt.close()
#     sp_receive_df.to_csv(out_file+'_receive_event.txt',sep='\t')
#     sp_donor_df.to_csv(out_file+'_donor_event.txt',sep='\t')

if 'calculate based on gene':
    sp_receive_gene_dic, sp_donor_gene_dic = {}, {}
    with open(out_file, 'r') as of:
        for i in of:
            cur_line=i.strip().split('\t')
            '''
            out file format example:
            donor_clade, receive_clade, donor_sp, receive_sp
            GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
            '''
            gf_id, donor_clade, receive_clade, donor_sp, donor_gene, receive_sp, receive_gene = cur_line[0], cur_line[4],cur_line[5],cur_line[6],cur_line[7],cur_line[8],cur_line[9]
            same_gene_trans = []
            for gene in donor_gene.split('|'):
                sp = gene.split('_')[0]
                gene_name_like = gene.rsplit('_', 1)[0]
                if gene_name_like in same_gene_trans and 'Lmul' not in gene_name_like and \
                        'Svir' not in gene_name_like and 'Sita' not in gene_name_like and \
                        'Bdis' not in gene_name_like and 'Pten' not in gene_name_like and \
                        'Asem' not in gene_name_like and 'Amyo' not in gene_name_like and \
                        'Dlat' not in gene_name_like and 'Amyo' not in gene_name_like and \
                        'Etef' not in gene_name_like and 'Sspo' not in gene_name_like:
                    continue
                same_gene_trans.append(gene_name_like)
                if sp in sp_receive_gene_dic:
                    for clade in receive_clade.split('|'):
                        if clade in sp_receive_gene_dic[sp]:
                            sp_receive_gene_dic[sp][clade] += 1
                        else:
                            sp_receive_gene_dic[sp][clade] = 1
                else:
                    for clade in receive_clade.split('|'):
                        if sp in sp_receive_gene_dic:
                            sp_receive_gene_dic[sp][clade] = 1
                        else:
                            sp_receive_gene_dic[sp] = {clade:1}
            same_gene_trans = []
            for gene in receive_gene.split('|'):
                sp = gene.split('_')[0]
                gene_name_like = gene.rsplit('_', 1)[0]
                if gene_name_like in same_gene_trans and 'Lmul' not in gene_name_like and \
                        'Svir' not in gene_name_like and 'Sita' not in gene_name_like and \
                        'Bdis' not in gene_name_like and 'Pten' not in gene_name_like and \
                        'Asem' not in gene_name_like and 'Amyo' not in gene_name_like and \
                        'Dlat' not in gene_name_like and 'Amyo' not in gene_name_like and \
                        'Etef' not in gene_name_like and 'Sspo' not in gene_name_like:
                    continue
                same_gene_trans.append(gene_name_like)
                if sp in sp_donor_gene_dic:
                    for clade in donor_clade.split('|'):
                        if clade in sp_donor_gene_dic[sp]:
                            sp_donor_gene_dic[sp][clade] += 1
                        else:
                            sp_donor_gene_dic[sp][clade] = 1
                else:
                    for clade in donor_clade.split('|'):
                        if sp in sp_donor_gene_dic:
                            sp_donor_gene_dic[sp][clade] = 1
                        else:
                            sp_donor_gene_dic[sp] = {clade:1}

    sp_receive_df = pd.DataFrame.from_dict(sp_donor_gene_dic, orient='index')
    sp_donor_df = pd.DataFrame.from_dict(sp_receive_gene_dic, orient='index')
    sp_receive_df['Clade'] = sp_receive_df.index.map(subg_clade_dic)
    sp_donor_df['Clade'] = sp_donor_df.index.map(subg_clade_dic)
    sp_receive_df = sp_receive_df.reindex(sp_ord)
    sp_donor_df = sp_donor_df.reindex(sp_ord)
    sp_receive_df = sp_receive_df.fillna(0)
    sp_donor_df = sp_donor_df.fillna(0)

    plt.figure(figsize=(6, 16))
    sns.set_style('whitegrid')

    group_gap = 1.5
    y_positions = []
    current_position = 0
    for group in sp_receive_df['Clade'].unique():
        group_data = sp_receive_df[sp_receive_df['Clade'] == group]
        y_positions.extend(np.arange(current_position, current_position + len(group_data)))
        current_position += len(group_data) + group_gap

    last_value = [0] * len(sp_receive_df)
    for group in sp_receive_df['Clade'].unique():
        plt.barh(y_positions, sp_receive_df[group].to_list(), left=last_value, linewidth=0, 
                color = template_color[group], edgecolor='none')
        last_value = [x+last_value[idx] for idx,x in enumerate(sp_receive_df[group].to_list())]

    last_value = [0] * len(sp_receive_df)
    for group in sp_donor_df['Clade'].unique():
        plt.barh(y_positions, [-num for num in sp_donor_df[group].to_list()], left=last_value, linewidth=0, 
                color=template_color[group], label=group, edgecolor='none')
        last_value = [-num+last_value[idx] for idx,num in enumerate(sp_donor_df[group].to_list())]


    xmin_value = sp_donor_df[list(template_color.keys())].sum(axis=1).max()
    xmax_value = sp_receive_df[list(template_color.keys())].sum(axis=1).max()
    ax = plt.gca()
    current_position = -0.5
    anno_rectangle = 50
    for group in sp_receive_df['Clade'].unique():
        group_data = sp_receive_df[sp_receive_df['Clade'] == group]
        width = len(group_data)
        rect = patches.Rectangle((-xmin_value-anno_rectangle, current_position), anno_rectangle, width,
                                linewidth=0.5, edgecolor='black', facecolor=template_color[group], alpha=0.6)
        ax.add_patch(rect)
        ax.text(-xmin_value-anno_rectangle/2, current_position+width/2, group,
                rotation=-90, horizontalalignment='center', verticalalignment='center', fontweight='bold')
        current_position += len(group_data) + group_gap

    plt.yticks(y_positions, list(sp_receive_df.index.values))
    # plt.ylabel('Value', labelpad=15)
    ax.set_xlim(-xmin_value-20, xmax_value)  # 限制 x 轴范围
    ax.set_ylim(-1, 128)
    plt.axvline(0, color='white', linewidth=1)
    plt.legend()
    plt.gca().xaxis.set_major_locator(MaxNLocator(10))
    sns.despine(left=True, bottom=True)
    plt.legend(loc='upper right')
    plt.tight_layout()
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\HGT_newcollection\3_outfile_stat\AGFD_gene_stat.pdf', bbox_inches='tight', format='pdf')
    plt.close()