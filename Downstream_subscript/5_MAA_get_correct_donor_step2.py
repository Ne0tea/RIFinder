'''
Descripttion: 根据Step1结果获得其中最小kaks值基因，制作id和homolog文件
Author: Ne0tea
version: 
Date: 2023-11-02 15:16:43
LastEditors: Ne0tea
LastEditTime: 2023-11-03 15:28:59
'''

import sys
import os
from random import sample
import re
import pandas as pd

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

def main(trans_file:str):
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
        give_gene={}
        donor_gene={}
        donor_kaks={}
        gene_99=[]
        result={}
        with open('./result/'+og_id+"_"+str(loc)+'/'+og_id+"_"+str(loc)+'_para/all-results.txt','r') as kfile:
            for i in kfile:
                if i.startswith('Seq'):
                    continue
                line=i.strip().split(' ')
                recive_gene=line[0].split('-')[0]
                gene=line[0].split('-')[1]
                ka=float(line[1])
                ks=float(line[2])
                kaks=float(line[3])
                if define_sub(recive_gene)==define_sub(gene) or define_sub(gene)=='OUT':
                    continue
                if gene in give_gene:
                    give_gene[gene].append(kaks)
                else:
                    give_gene[gene]=[kaks]
                if ka > 20 or ks>20 :
                    gene_99.append(gene)
        for i in give_gene:
            if i in gene_99:
                continue
            else:
                result[i]=sum(give_gene[i])/len(give_gene[i])
        min_gene=min(zip(result.values(), result.keys()))
        print(HGT_gene[0],min_gene[1],min_gene[0])
        with open('./gene_id_and_homo/'+og_id+"_"+str(loc)+"ddd.homolog",'w') as homo_file:
            for x in HGT_gene:
                homo_file.write(str(min_gene[1])+"\t"+str(x)+"\n")
            # homo_file.write(str(out_gene)+"\t"+str(donor_gene)+"\n")
        with open('./gene_id_and_homo/'+og_id+"_"+str(loc)+"ddd_gene_id",'w') as gene_file:
            for x in HGT_gene:
                gene_file.write(str(x)+"\n")
            gene_file.write(str(min_gene[1])+"\n")



if __name__ == "__main__":
    # trans_file=sys.argv[1]
    trans_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\0724_manully_check.xlsx"
    main(trans_file)