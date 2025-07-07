'''
Descripttion: 循环计算HGT中一个gene对应在树中的所有kaks值
Author: Ne0tea
version: 
Date: 2023-11-02 15:16:43
LastEditors: Ne0tea
LastEditTime: 2023-11-03 15:29:00
'''

from ete3 import Tree
import os
from random import sample
import re
import pandas as pd

def main(trans_file:str,tree_fp:str):
    if not os.path.isdir('./gene_id_and_homo'):
        os.mkdir('./gene_id_and_homo')
    trans=pd.read_excel(trans_file,engine='openpyxl',sheet_name=1)
    pattern = r"'([^'\s]*)'"
    for loc,row in trans.iterrows():
        # print(loc)
        if row['Confirm']:
            signal=row['Type']
            og_id=row['OG_ID']
            Gene_list=row['Gene']
            if signal=="G_line_R":
                HGT_gene=Gene_list.split("|")
            elif signal=="G_line_A":
                HGT_gene=re.findall(pattern,Gene_list)
                # HGT_gene="|".join(HGT_gene)
        else:
            continue
        give_gene=[]
        tree_file=os.path.join(tree_fp,og_id+"_ali.mafft_rmdup_trimal_JTT.treefile")
        gene_tree:Tree=Tree(tree_file)
        all_gene=gene_tree.get_leaf_names()
        sus_donor_gene=list(set(all_gene) - set(HGT_gene))
        target_gene=sample(HGT_gene, min(3,len(HGT_gene)))
        with open('./gene_id_and_homo/'+og_id+"_"+str(loc)+".homolog",'w') as homo_file:
            for i in target_gene:
                for x in sus_donor_gene:
                    homo_file.write(str(i)+"\t"+str(x)+"\n")
            # homo_file.write(str(out_gene)+"\t"+str(donor_gene)+"\n")
        with open('./gene_id_and_homo/'+og_id+"_"+str(loc)+"_gene_id",'w') as gene_file:
            for x in sus_donor_gene:
                gene_file.write(str(x)+"\n")
            for x in target_gene:
                gene_file.write(str(x)+"\n")



if __name__ == "__main__":
    # trans_file=sys.argv[1]
    #tree_fp=sys.argv[2]
    tree_fp=r'E:\Bio_analysis\Grass_horizatal_transversion\7_final_plot\0805manully_check_circular_tree_set'
    trans_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\0724_manully_check.xlsx"
    main(trans_file,tree_fp)