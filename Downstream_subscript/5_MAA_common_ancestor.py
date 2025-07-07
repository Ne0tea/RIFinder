'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-10-18 10:58:29
LastEditors: Ne0tea
LastEditTime: 2023-10-18 11:43:29
'''
from ete3 import Tree
import pandas as pd
import re


def main(trans_file:str,tree:str):
    sp_tree=Tree(tree)
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
        # print("|".join(HGT_gene),end='\t')
        print(og_id,loc,sep='_',end='\t')
        print(len(HGT_gene),end='\t')
        HGT_sp=[x.split('_')[0] for x in HGT_gene]
        if len(set(HGT_sp))==1:
            common_ancestor_sp=list(set(HGT_sp))
        else:
            common_ancestor_sp=sp_tree.get_common_ancestor(list(set(HGT_sp)))
        print(len(common_ancestor_sp))

if __name__ == "__main__":
    # trans_file=sys.argv[1]
    #sp_tree=sys.argv[2]
    trans_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\0724_manully_check.xlsx"
    sp_tree=r"E:\Bio_analysis\Grass_horizatal_transversion\poaceae_sptree_rooted.nwk"
    main(trans_file,sp_tree)