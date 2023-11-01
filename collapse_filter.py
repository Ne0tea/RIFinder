'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-06-02 13:47:26
LastEditors: Ne0tea
LastEditTime: 2023-10-20 23:14:18
'''
import sys
import logging
from ete3 import Tree
import global_variable as glv

def define_sub(gene):
    clade_subg_dic=glv.get('clade_subg_dic')
    sp=gene.split("_")[0]
    for i in clade_subg_dic:
        if sp in clade_subg_dic[i]:
            return i

def collapsed_leaf(node):
    sub_list=[]
    for i in node2labels1[node]:
        # print(i.name)
        if i != "":
            sub=define_sub(i)
            if not sub :
                logging.warning('Can\'t find according sub clade for {0}'.format(i))
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

def main(tree_file:str):
    new_tree=Tree(tree_file,format=1)
    global node2labels1
    node2labels1 = new_tree.get_cached_content(store_attr="name")
    collapse_tree = Tree(new_tree.write(is_leaf_fn=collapsed_leaf))
    sub_dic={}
    sub_suc={}
    sub_stat={}
    for clade in collapse_tree:
        up=clade.up
        if up.get_sisters() and clade.name !="OUT":
            if clade.name not in sub_dic:
                sub_dic[clade.name]=1
                if clade.name == up.get_sisters()[0].name:
                    # sub_dic[clade.name]+=1
                    if clade.name not in sub_suc:
                        sub_suc[clade.name]=1
                    else:
                        sub_suc[clade.name]+=1
                else:
                    sub_dic[clade.name]+=1
            else:
                if clade.name == up.get_sisters()[0].name:
                    sub_dic[clade.name]+=1
                    if clade.name not in sub_suc:
                        sub_suc[clade.name]=1
                    else:
                        sub_suc[clade.name]+=1
                else:
                    sub_dic[clade.name]+=1

    # og=re.search('OG[0-9]{7}',tree_file).group()
    # print(og,end="\t")
    statue=True
    count=0
    reason=''
    for i in sorted(sub_dic.keys()):
        if i not in sub_suc:
            i_suc=0
        else:
            i_suc=sub_suc[i]
        # sub_stat[i]=i_suc/sub_dic[i]
        if i_suc/sub_dic[i] < 0.3:
            count+=sub_dic[i]
        if count > 8:
            statue=False
        reason=reason+str(i_suc/sub_dic[i])+"/"+str(sub_dic[i])
    return reason,statue

if __name__ == "__main__":
    # tree_file=sys.argv[1]
    # tf_file=sys.argv[2]
    tree_file=r"D:\study resource\HGT\4_significant_gene\OG0003936_ali.mafft_rmdup_trimal_JTT.treefile"
    # tf_file=r"D:\Bio_analysis\Grass_horizatal_transversion\test\costume_check\OG0000568_transfer_pair"
    reason,state=main(tree_file)
    print(reason,state)


# BB=["Phedu","Dlat","Olat","Bamp","DlatA","DlatB","DlatC"]
# #OO=["Lper","Oruf","Osat","Obra","Omer","Oglu","Opun","Ogla","Obar","Zlat","Zpal"]
# OO=["Lper","Oruf","Osat","Obra","Omer","Oglu","Opun","Ogla","Obar","Zlat","Zpal"]
# #PP=["Aspl","Atau","Astr","HvulSp","Pten","Scer","Telo","Ttur","Bdis","Bsta","HvulC","Taes","Tdic","Tura"]
# PP=["Aspl","Atau","Astr",'Hvul', "HvulSp","Pten","Scer","Telo","Ttur","TturA","TturB","Bdis","Bsta","HvulC",\
#   "Taes","TaesA","TaesB","TaesD","Tdic","TdicA","TdicB","Tura"]
# #PAC=["Ehap","Asem","Came","Cpur","Doli","Dexi","Msin","Phal","Pvir","Sita","Svir","Sbic","Zmay","Ecol","Ecru","Eory"]
# PAC=["Ehap","Asem","Came","Cpur","CpurA","CpurB","Doli","Dexi","DexiA","DexiB","Msin","MsinA","MsinB","Phal",\
#   "Pvir","PvirK","PvirN","Sita","Svir","Sbic","Zmay",'EcolDL',"Ecol","EcolD","EcolFL","EcolE","EcolEL","EcolF","Ecru",\
#     "EcruBH","EcruA","EcruB","EcruC","Eory","EoryA","EoryB"]
# #MAD=["Cson","Ecur","Etef","Enin","Otho","Zjap"]
# MAD=["Cson","CsonA","CsonB","Ecur","EtefUn","Etef","EtefA","EtefB","Enin","Otho","Zjap"]