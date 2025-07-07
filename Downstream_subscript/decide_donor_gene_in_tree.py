'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-08-16 15:33:05
LastEditors: Ne0tea
LastEditTime: 2023-10-28 20:43:28
'''
from ete3 import Tree
import os
import statistics
import filter_length
import sys
import re
import pandas as pd
from collections import Counter
###亚科列表
BB=["Phedu","Dlat","Olat","Bamp","DlatA","DlatB","DlatC"]
#OO=["Lper","Oruf","Osat","Obra","Omer","Oglu","Opun","Ogla","Obar","Zlat","Zpal"]
OO=["Lper","Oruf","Osat","Obra","Omer","Oglu","Opun","Ogla","Obar","Zlat","Zpal"]
#PP=["Aspl","Atau","Astr","HvulSp","Pten","Scer","Telo","Ttur","Bdis","Bsta","HvulC","Taes","Tdic","Tura"]
PP=["Aspl","Atau","Astr",'Hvul', "HvulSp","Pten","Scer","Telo","Ttur","TturA","TturB","Bdis","Bsta","HvulC",\
  "Taes","TaesA","TaesB","TaesD","Tdic","TdicA","TdicB","Tura"]
#PAC=["Ehap","Asem","Came","Cpur","Doli","Dexi","Msin","Phal","Pvir","Sita","Svir","Sbic","Zmay","Ecol","Ecru","Eory"]
PAC=["Ehap","Asem","Came","Cpur","CpurA","CpurB","Doli","Dexi","DexiA","DexiB","Msin","MsinA","MsinB","Phal",\
  "Pvir","PvirK","PvirN","Sita","Svir","Sbic","Zmay",'EcolDL',"Ecol","EcolD","EcolFL","EcolE","EcolEL","EcolF","Ecru",\
    "EcruBH","EcruA","EcruB","EcruC","Eory","EoryA","EoryB"]
#MAD=["Cson","Ecur","Etef","Enin","Otho","Zjap"]
MAD=["Cson","CsonA","CsonB","Ecur","EtefUn","Etef","EtefA","EtefB","Enin","Otho","Zjap"]
OUT=["Plat","Acom"]
sp_list=["Phed","Dlat","Olat","Bamp","Dlat","Lper","Oruf","Osat","Obra","Omer","Oglu","Opun","Ogla","Obar","Zlat","Zpal",\
    "Aspl","Atau","Astr",'Hvul',"Pten","Scer","Telo","Ttur","Bdis","Bsta","Eory","Cson",\
  "Taes","Tdic","Tura","Ehap","Asem","Came","Cpur","Doli","Dexi","Msin","Phal","Enin","Otho",\
  "Pvir","Sita","Svir","Sbic","Zmay","Ecol","Ecru","Ecur","Etef","Zjap"]
###亚科颜色列表
color_template={"BB":"#b71515","OO":"#e97a0c","PP":"#ffde0a","PAC":"#034732","MAD":"#092e86","OUT":"#a9a29c"}
tansfer_tmplate={"BB":"#b71515","OO":"#e97a0c","PP":"#ffde0a","PAC":"#034732","MAD":"#092e86"}

def define_sub(gene):
    sp=gene.split("_")[0]
    if sp in BB:
        return "BB"
    elif sp in OO:
        return "OO"
    elif sp in PP:
        return "PP"
    elif sp in PAC:
        return "PAC"
    elif sp in MAD:
        return "MAD"
    elif sp in OUT:
        return "OUT"
    else:
        return "empty"

def get_give_gene(gene:str,gene_tree:Tree):
    give_list=[]
    cur_node=gene_tree&gene
    cur_sub=define_sub(cur_node.name)
    try:
        check_node=cur_node.up.get_sisters()[0]
    except IndexError:
        return
    pos_give=[define_sub(i) for i in check_node.get_leaf_names()]
    pos_give_gene=[i for i in check_node.get_leaf_names()]
    while len(set(pos_give))==1 and list(set(pos_give))[0]==cur_sub:
        check_node=check_node.up
        check_node=check_node.get_sisters()[0]
        pos_give=[define_sub(i) for i in check_node.get_leaf_names()]
        pos_give_gene=[i for i in check_node.get_leaf_names()]
    if len(pos_give) > 4:
        if len(set(pos_give))==1:
            give_sub=pos_give[0]
            give_list.extend(pos_give_gene)
        else:
            result=Counter(pos_give)
            # print(result)
            give_sub=max(result.keys(), key=result.get)
            # print(give_sub)
            if give_sub==cur_sub:
                result.pop(give_sub)
                while max(result.values())<int(4):
                    check_node=check_node.up.get_sisters()[0]
                    pos_give=[define_sub(i) for i in check_node.get_leaf_names() if define_sub(i)!=cur_sub]
                    result.update(pos_give)
                    pos_give_gene.extend([i for i in check_node.get_leaf_names() if define_sub(i)!=cur_sub])
                give_sub=max(result.keys(), key=result.get)
            give_list.extend([i for i in pos_give_gene if define_sub(i)==give_sub])
    else:
        result=Counter(pos_give)
        while max(result.values())<int(4) and check_node.up.get_sisters():
            # print(gene,result)
            check_node=check_node.up.get_sisters()[0]
            pos_give=[define_sub(i) for i in check_node.get_leaf_names() if define_sub(i)!=cur_sub]
            result.update(pos_give)
            pos_give_gene.extend([i for i in check_node.get_leaf_names() if define_sub(i)!=cur_sub])
        give_sub=max(result.keys(), key=result.get)
        give_list.extend([i for i in pos_give_gene if define_sub(i)==give_sub])
    return give_list


def main(trans_file:str,tree_fp:str):
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
        # if og_id !='OG0000776':
        #     continue
        give_gene=[]
        tree_file=os.path.join(tree_fp,og_id+"_ali.mafft_rmdup_trimal_JTT.treefile")
        gene_tree:Tree=Tree(tree_file)
        outleaf=filter_length.find_outlength(gene_tree)
        intree_gene=[]
        for i in HGT_gene :
            # print(i)
            try:
                cur_node=gene_tree&i
            except :
                print("there is no",i,"in",og_id)
                continue
            intree_gene.append(i)
            cur_sub=define_sub(cur_node.name)
            if cur_sub=="OUT":
                continue
            cur_set=get_give_gene(i,gene_tree)
            if cur_set is not None:
                give_gene.extend(cur_set)
                give_gene=list(set(give_gene))
            else:
                continue
        give_gene_dist=dict((x,statistics.median([(gene_tree&x).get_distance(i) for i in intree_gene]))for x in give_gene )
        donor_gene=min(give_gene_dist.keys(),key=give_gene_dist.get)
        out_gene_dic=dict((i,(gene_tree&donor_gene).get_distance(i)) for i in gene_tree.get_leaf_names() if i not in HGT_gene and i != donor_gene)
        out_gene=min(out_gene_dic.keys(),key=out_gene_dic.get)
        with open('xx','a') as result:
            line='\t'.join([og_id+'_'+str(loc),donor_gene,out_gene,define_sub(donor_gene),cur_sub])
            # print(og_id+'_'+str(loc),donor_gene,out_gene,define_sub(donor_gene),cur_sub,sep='\t')
            result.write(line+'\n')

        # with open(og_id+"_"+str(loc)+".homolog",'w') as homo_file:
        #     for i in HGT_gene:
        #         homo_file.write(str(out_gene)+"\t"+str(i)+"\n")
        #     homo_file.write(str(out_gene)+"\t"+str(donor_gene)+"\n")
        # with open(og_id+"_"+str(loc)+"_gene_id",'w') as gene_file:
        #     for i in HGT_gene:
        #         gene_file.write(str(i)+"\n")
        #     gene_file.write(str(donor_gene)+"\n")
        #     gene_file.write(str(out_gene)+"\n")


if __name__ == "__main__":
    # trans_file=sys.argv[1]
    #tree_fp=sys.argv[2]
    tree_fp=r'E:\Bio_analysis\Grass_horizatal_transversion\7_final_plot\0805manully_check_circular_tree_set'
    trans_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\0724_manully_check.xlsx"
    main(trans_file,tree_fp)