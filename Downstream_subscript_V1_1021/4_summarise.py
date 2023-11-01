'''
Descripttion: INPUT: 鉴定HGT结果 OUTPUT: csv结尾矩阵文件,需要转换为xlsx实现可视化(复制粘贴即可) 未删除HGT中长枝基因
Author: Ne0tea
version:v2
Date: 2023-06-13 19:07:28
LastEditors: Ne0tea
LastEditTime: 2023-10-24 13:52:49
'''
import re
import sys
import os
import pandas as pd
import numpy as np
from ete3 import Tree
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
sp_list=["Phedu","Dlat","Olat","Bamp","Dlat","Lper","Oruf","Osat","Obra","Omer","Oglu","Opun","Ogla","Obar","Zlat","Zpal",\
    "Aspl","Atau","Astr",'HvulC','HvulSp',"Pten","Scer","Telo","Ttur","Bdis","Bsta","Eory","Cson",\
  "Taes","Tdic","Tura","Ehap","Asem","Came","Cpur","Doli","Dexi","Msin","Phal","Enin","Otho",\
  "Pvir","Sita","Svir","Sbic","Zmay","Ecol","Ecru","Ecur","Etef","Zjap"]
###亚科颜色列表
color_template={"BB":"#b71515","OO":"#e97a0c","PP":"#ffde0a","PAC":"#034732","MAD":"#092e86","OUT":"#a9a29c"}
tansfer_tmplate={"BB":"#b71515","OO":"#e97a0c","PP":"#ffde0a","PAC":"#034732","MAD":"#092e86"}

gene_len={}
bed_file=r'E:\Bio_analysis\Grass_horizatal_transversion\ful_gff.bed'
with open (bed_file,'r') as bedf:
    for i in bedf.readlines():
        c_line=i.strip().split()
        gene_id=c_line[1]
        gene_s=c_line[2]
        gene_e=c_line[3]
        lennn=int(gene_e)-int(gene_s)
        gene_len[gene_id]=lennn

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

def get_give_sub(gene:str,gene_tree:Tree):
    sp=gene[0:4]
    # print(sp)
    homo_list=[x for x in gene_tree.get_leaf_names() if x[0:4]==sp]
    cur_node=gene_tree&gene
    cur_sub=define_sub(cur_node.name)
    check_node=cur_node.up.get_sisters()[0]
    pos_give=[define_sub(i) for i in check_node.get_leaf_names()]
    while len(set(pos_give))==1 and list(set(pos_give))[0]==cur_sub:
        check_node=check_node.up
        check_node=check_node.get_sisters()[0]
        pos_give=[define_sub(i) for i in check_node.get_leaf_names()]
    if len(pos_give) > 4:
        if len(set(pos_give))==1:
            give_sub=pos_give[0]
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
                give_sub=max(result.keys(), key=result.get)
    else:
        result=Counter(pos_give)
        while max(result.values())<int(4) and check_node.up.get_sisters():
            # print(result)
            check_node=check_node.up.get_sisters()[0]
            pos_give=[define_sub(i) for i in check_node.get_leaf_names() if define_sub(i)!=cur_sub]
            result.update(pos_give)
        give_sub=max(result.keys(), key=result.get)
    return give_sub

def get_give_sub2(gene:str,gene_tree:Tree):
    cur_node=gene_tree&gene
    cur_sub=define_sub(cur_node.name)
    check_node=cur_node.up.get_sisters()[0]
    check_node=[i for i in check_node if i.name != '' ][0]
    give_sub=[]
    switch=0
    while switch < 2:
        check_sub=define_sub(check_node.name)
        # print(check_sub)
        if check_sub=="empty":
            children_list=check_node.get_leaves()
            children_sub=[define_sub(i.name) for i in children_list]
            result=Counter(children_sub)
            increase=1
            for x in result:
                if x != cur_sub and increase < result[x]:
                    switch-=increase
                    increase=result[x]
                    switch+=increase
                    give_sub.append(x)
                    # print(i,sp,cur_sub,x)
                    break
        else:
            if check_sub!=cur_sub:
                switch+=1
                give_sub.append(check_sub)
                # print(i,sp,cur_sub,check_sub)

        check_node=check_node.up
        check_node=check_node.get_sisters()[0]
    # print(len(set(give_sub)))
    if len(set(give_sub))==1:
        give_sub=give_sub[0]
    else:
        give_sub=give_sub[1]
    if give_sub == "OUT":
        return
    return give_sub

def main(trans_file:str,gene_file_path:str,outpt:dir,manul_check=False,sheet=0):
    # items = os.listdir(gene_file_path)
    # tree_file=r"D:\Bio_analysis\Grass_horizatal_transversion\3_KNN\OG0000856_ali.mafft_rmdup_trimal_JTT.treefile"
    # trans_file=r'D:\study resource\HGT\0721_gene'
    global gene_len
    summary_out=os.path.join(outpt, 'transfer_summary.txt')
    matrix_out=os.path.join(outpt, 'transfer_matrix.csv')
    # df = pd.read_excel(trans_file,engine='openpyxl')
    df_matrix=pd.DataFrame(columns=[i+"_2_"+y for i in tansfer_tmplate for y in tansfer_tmplate if y != i],index=[sp_list])
    trans_dic={}
    pattern = r"'(.*?)'"

#获得OG为键的字典
    if manul_check:
        df = pd.read_excel(trans_file,engine='openpyxl',sheet_name=sheet)
        len_set=[]
        avg_set=[]
        count_set=[]
        for _,row in df.iterrows():
            tt_length=0
            if row['Confirm']:
                signal=row['Type']
                og_id=row['OG_ID']
                Gene_list=row['Gene']
                if signal=="G_line_R":
                    HGT_gene=Gene_list
                elif signal=="G_line_A":
                    HGT_gene=re.findall(pattern,Gene_list)
                    HGT_gene="|".join(HGT_gene)
                if og_id not in trans_dic:
                    trans_dic[og_id]=HGT_gene
                else:
                    trans_dic[og_id]=trans_dic[og_id]+'|'+HGT_gene
                for x in HGT_gene.split("|"):
                    tt_length+=gene_len[x]
                len_set.append(tt_length)
                avg_set.append(tt_length/len(HGT_gene.split("|")))
                count_set.append(len(HGT_gene.split("|")))
            else:
                len_set.append(0)
                avg_set.append(0)
        df['Len']=len_set
        df['avg']=avg_set
        df['count']=count_set
    else:
        with open(trans_file,'r') as tr_f:
            for i in tr_f.readlines():
                c_line=i.strip().split()
                signal=c_line[0]
                og_id=c_line[1]
                if signal=="G_line_R":
                    HGT_gene=c_line[2]
                elif signal=="G_line_A":
                    HGT_gene=re.findall(pattern,c_line[2])
                    HGT_gene="|".join(HGT_gene)
                if og_id not in trans_dic:
                    trans_dic[og_id]=HGT_gene
                else:
                    trans_dic[og_id]=trans_dic[og_id]+'|'+HGT_gene

#制作summary文件
    summary=open(summary_out,'w')
    # count=0
    for i in trans_dic:
        tree_file=os.path.join(gene_file_path,i+"_ali.mafft_rmdup_trimal_JTT.treefile")
        # print(tree_file)
        og=i
        transfer_set=trans_dic[og]
        try:
            tf_gene_list=transfer_set.split("|")
        except AttributeError:
            print("there is no ",og)
            continue
        tf_gene_list=[s.strip() for s in tf_gene_list]
        tf_gene_list=list(filter(None, tf_gene_list))
        # print(tf_gene_list)
        gene_tree=Tree(tree_file)
        for x in gene_tree:
            if define_sub(x.name)=="OUT":
                gene_tree.set_outgroup(x)
                break
        for y in tf_gene_list:
            if y.startswith('Hvul'):
                sp=y.split('_')[0]
            elif y.startswith('Phedu'):
                sp='Phedu'
            else:
                sp=y[0:4]
            try:
                cur_node=gene_tree&y
            except :
                print("there is no",y,"in",og)
                continue
            cur_sub=define_sub(cur_node.name)
            if cur_sub=="OUT":
                continue
            give_sub=get_give_sub(y,gene_tree)
            # give_sub2=get_give_sub2(y,gene_tree)
            # if give_sub2 != give_sub:
            #     print(og,y,sp,cur_sub,give_sub,give_sub2)
            line1=" ".join([og,y,sp,cur_sub,give_sub+"\n"])
            if df_matrix.loc[sp,give_sub+"_2_"+cur_sub] is np.nan:
                df_matrix.loc[sp,give_sub+"_2_"+cur_sub]=1
            else:
                # print(df_matrix.loc[sp,give_sub+"_2_"+cur_sub])
                df_matrix.loc[sp,give_sub+"_2_"+cur_sub]+=1
            # print(df_matrix)
            summary.write(line1)
            # summary.write(line2)

    summary.close()
    df_matrix.to_csv(matrix_out,sep=',',index=True,header=True)
    df.to_csv(trans_file.replace('.xlsx','_len.csv'),sep=',',index=True,header=True)

if __name__ == "__main__":
    # trans_file=sys.argv[1]
    # tree_fp=sys.argv[2]
    # workpt=sys.argv[3]
    # sheet=sys.argv[5]
    tree_fp=r'E:\Bio_analysis\Grass_horizatal_transversion\7_final_plot\0805manully_check_circular_tree_set'
    workpt=r'E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration'
    trans_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\0724_manully_check.xlsx"
    sheet=1
    if trans_file.endswith('xlsx'):
        manul_check=True
        # print(globals())
        if 'sheet' not in globals():
            print('no sheet assigned, it will calcualte based on sheet 0')
            sheet=0
    main(trans_file,tree_fp,workpt,manul_check,1)