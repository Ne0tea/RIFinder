'''
Descripttion: 鉴定在基因树中的clade转移,完全忽视single gene的HGT;
              INPUT:gene_tree ML NWK格式的基因树文件 
                    sup_value 过滤bootstrap值
              OUTPUT:print输出OG_id state HGT_gene_list(sep="\t") number_of_HGT_gene
Update:
    v2将树结构识别为区域,即为PACMAD块和OOBBPP块,检测两者之间的交叉
    V3增加了对单系的识别,具体例子为OG0000766的OO类群;
    V4增加了对单系是否识别为区域的判断;
    V5增加Bootstrap筛选;
    V6更新发生HGT结构判断标准,置于两亚科之间的HGT认定为Suspetic类型;
Author: Ne0tea
version: v7
Date: 2023-07-08 20:46:51
LastEditors: Ne0tea
LastEditTime: 2023-07-19 22:35:27
'''

import sys
from ete3 import Tree
import re
import logging
from filter_length import find_outlength
import global_variable as glv

###V4更新：如果clade中存在一个单系与上下的clade类型都不一致，则该clade需要进行分类
def define_HGT_clade(gene_list:list,gene_tree:Tree,calde:dict,standard:int=50,test_mode:bool=False) -> bool:
    
    sub_list=[i.split("_")[0] for i in gene_list]
    clade_list=[ list(calde.keys())[0] if i in calde[list(calde.keys())[0]]  else list(calde.keys())[1] for i in sub_list]
    if clade_list[0]==clade_list[1]:
        # print(11,clade_list)
        for i in range(2,len(clade_list)-1):
            support=(gene_tree&gene_list[i]).up.support
            # print(gene_list[i],support)
            if clade_list[i] != clade_list[i-1] and clade_list[i] != clade_list[i+1] and support > standard:#V5增加了根据Bootstrap的筛选
                # if sub_list[i] != sub_list[i-1] and sub_list[i] != sub_list[i+1] and sub_list[i-1] == sub_list[i+1]:#V5中进行了删除
                # print(support)
                if sub_list[i-1] == sub_list[i+1]:
                    if test_mode:
                        print(clade_list,'1 = 2')
                    return True, min(clade_list,key=clade_list.count)
                else:
                    if i==2 and sub_list[3] in [sub_list[0],sub_list[1]]:
                        return True, min(clade_list,key=clade_list.count)
                    else:
                        return "suspicious",min(clade_list,key=clade_list.count) #代表clade为可能发生HGT的clade(多数情况为PACMAD被包裹在BPO中，但是分类整齐)，存在clade分类的必要
        return False,None #列表中不存在与前后clade都不相同的leaf，说明不发生HGT
    else:
        support=(gene_tree&gene_list[0]).up.support
        if support > standard:
            # print(gene_list,support)
            if test_mode:
                print(clade_list,'1 != 2 and support > standard')
                print(sub_list)
            if clade_list[2] == list(calde.keys())[0]:
                return True,None
            elif sub_list[2] in [sub_list[0],sub_list[1]]:
                if test_mode:
                    print(sub_list,'sub 3 in sub1,2')
                return True,None
            else:
                return False,None
        else:
            for i in range(2,len(clade_list)-1):
                support=(gene_tree&gene_list[i]).up.support
                # print(gene_list[i],support)
                if clade_list[i] != clade_list[i-1] and clade_list[i] != clade_list[i+1] and support > standard:#V5增加了根据Bootstrap的筛选
                    # if sub_list[i] != sub_list[i-1] and sub_list[i] != sub_list[i+1] and sub_list[i-1] == sub_list[i+1]:#V5中进行了删除
                    # print(support)
                    if sub_list[i-1] == sub_list[i+1]:
                        if test_mode:
                            print(clade_list,'1 != 2 and support < standard')
                            print(sub_list)
                        return True, min(clade_list,key=clade_list.count)
                    else:
                        return "suspicious",min(clade_list,key=clade_list.count)
            return False,None #列表中不存在与前后clade都不相同的leaf，说明不发生HGT
###
'''
根据gene_name确定科属 
INPUT:gene_name:str
OUTPUT:sub_family:str
'''
def define_sub(gene):
    clade_subg_dic=glv.get('clade_subg_dic')
    sp=gene.split("_")[0]
    for i in clade_subg_dic:
        if sp in clade_subg_dic[i]:
            return i

'''
根据基因树中的外类群对基因树重定根
INPUT:gene_tree:TREE
OUT:TREE:list 重定根后的树对象列表
'''
def re_root(gene_tree:Tree) -> Tree:
    outsp=[i for i in gene_tree.get_leaves if define_sub(i) == "OUT" ]
    reroot_tree_list=[]
    for i in outsp:
        reroot_tree=gene_tree.set_outgroup(gene_tree&i)
        reroot_tree_list.append(reroot_tree)
    if not outsp:
        reroot_tree=gene_tree.get_midpoint_outgroup()
        reroot_tree_list.append(reroot_tree)
    return reroot_tree_list

'''
根据collapsed基因树,过滤长枝之后,鉴定古老HGT产生的与物种树拓扑矛盾的高bootstrap值的HGT基因集:
INPUT:collapsed_gene_tree:Tree
      standard:int
OUTPUT:print输出OG_id state HGT_gene_list(sep="\t") number_of_HGT_gene
'''

def find_clade(gene_tree:Tree,standard:int=50,test_mode:bool=False):

    sub_phlo=glv.get('sub_phlo')
    phylo=Tree(sub_phlo)
    calde={}
    if 'OUT' in phylo.get_leaf_names():
        main=[i for i in phylo.get_children() if i.get_leaf_names() != "OUT" ]
    else:
        main=phylo
    for i in main.get_children():
        calde["".join(i.get_leaf_names())]=i.get_leaf_names()
    HGT_list=[]
    sus_list=[]
    sister_list=[]
    double_check=[]
    for i in gene_tree:
        # print(i)
        try:
            sister=i.get_sisters()[0]
        except IndexError:
            continue
        if sister.name != "" and i.name not in sister_list:
            sister_list.append(sister.name)
            clade_list=[i.name,sister.name]
            try:
                up_sister=i.up.get_sisters()[0]
            except IndexError:
                continue
            while up_sister.name != "" :
                # print(up_sister.name)
                clade_list.append(up_sister.name)
                try:
                    up_sister=up_sister.up.get_sisters()[0]
                except IndexError:
                    break
            # print(up_sister)
            clade_dic={list(calde.keys())[0]:[],list(calde.keys())[1]:[]}
            PACMAD=0
            BPO=0
            if len(clade_list)==2:
                double_check.append(i)
                continue
            # print(clade_list)
            statue,clade_HGT=define_HGT_clade(clade_list,gene_tree,calde,standard=50,test_mode=test_mode)
            if test_mode:
                print(statue)
            if statue :
                for x in clade_list:
                    sub=x.split("_")[0]
                    if sub in calde[list(calde.keys())[0]]:
                        PACMAD+=int(x.split("_")[1])
                        clade_dic[list(calde.keys())[0]].append(x)
                    elif sub in calde[list(calde.keys())[1]]:
                        BPO+=int(x.split("_")[1])
                        clade_dic[list(calde.keys())[1]].append(x)
                # print(BPO,PACMAD,clade_HGT) 
                if PACMAD > BPO and (not clade_HGT or clade_HGT==list(calde.keys())[1]):#V3中clade的分类标准;V5中增加后半段
                    if statue != "suspicious":
                        HGT_list.extend(clade_dic[list(calde.keys())[1]])
                    else:
                        sus_list.extend(clade_dic[list(calde.keys())[1]])
                    up_sister.get_sisters()[0].add_features(sub=list(calde.keys())[0])
                    up_sister.get_sisters()[0].add_features(sub_list="_".join(list(set(clade_list)&set(list(calde.keys())[0]))))
                    # print(up_sister.get_sisters()[0])
                elif BPO > PACMAD and (not clade_HGT or clade_HGT==list(calde.keys())[0]):#V3中clade的分类标准;V5中增加后半段
                    if statue != "suspicious":
                        HGT_list.extend(clade_dic[list(calde.keys())[0]])
                    else:
                        sus_list.extend(clade_dic[list(calde.keys())[0]])    
                    up_sister.get_sisters()[0].add_features(sub=list(calde.keys())[1])
                    up_sister.get_sisters()[0].add_features(sub_list="_".join(list(set(clade_list)&set(list(calde.keys())[1]))))
            else:
                for x in clade_list:
                    sub=x.split("_")[0]
                    if sub in calde[list(calde.keys())[0]]:
                        PACMAD+=int(x.split("_")[1])
                        clade_dic[list(calde.keys())[0]].append(x)
                    elif sub in calde[list(calde.keys())[1]]:
                        BPO+=int(x.split("_")[1])
                        clade_dic[list(calde.keys())[1]].append(x)
                if BPO==0 :
                    up_sister.get_sisters()[0].add_features(sub=list(calde.keys())[0])
                    up_sister.get_sisters()[0].add_features(sub_list="_".join(list(set(clade_list)&set(calde[list(calde.keys())[0]]))))
                elif PACMAD==0:
                    up_sister.get_sisters()[0].add_features(sub=list(calde.keys())[1])
                    up_sister.get_sisters()[0].add_features(sub_list="_".join(list(set(clade_list)&set(calde[list(calde.keys())[1]]))))
                # print(up_sister.get_sisters()[0])
###V3更新
    # print(double_check)
    for i in double_check:
        # print(i.name)
        sister=i.get_sisters()[0]
        # print(sister.up.get_sisters()[0])
        sub=list(calde.keys())[0] if i.name.split("_")[0] in calde[list(calde.keys())[0]] else list(calde.keys())[1]
        sister_sub=list(calde.keys())[0] if sister.name.split("_")[0] in calde[list(calde.keys())[0]] else list(calde.keys())[1]
        if sub == sister_sub:
            i.up.add_features(sub=sub)
    for i in double_check:
        # print(i.name)
        sister=i.get_sisters()[0]
        # print(sister.up.get_sisters()[0])
        sub=list(calde.keys())[0] if i.name.split("_")[0] in calde[list(calde.keys())[0]] else list(calde.keys())[1]
        sister_sub=list(calde.keys())[0] if sister.name.split("_")[0] in calde[list(calde.keys())[0]] else list(calde.keys())[1]
        # print(i.up.get_sisters()[0],sub)
        try:
            # if i.up.get_sisters()[0].sub == sub and sub != sister_sub:
#V5增加了根据Bootstrap的筛选
            if i.up.get_sisters()[0].sub == sub and sub != sister_sub and i.up.support > standard \
                and i.name.split("_")[0] in i.up.get_sisters()[0].sub_list.split("_"):
                HGT_list.append(sister.name)
                if test_mode:
                    print(i,'clade len 2 and success')
            elif i.up.get_sisters()[0].sub == sister_sub and sub != sister_sub and i.up.support > standard \
                and sister.name.split("_")[0] in i.up.get_sisters()[0].sub_list.split("_"):
                if test_mode:
                    print(i,'clade len 2 and success')
                HGT_list.append(i.name)
        except AttributeError:
            continue
###
    return HGT_list,sus_list

'''
树结构,分布在相同拓扑结构中的sub_gene合并为同一枝
'''
def collapsed_leaf(node:Tree) -> bool:
    global pick
    sub_list=[]
    gene_list=[]
    # gene_dic[node.get_topology_id()]=[]
    for i in node2labels1[node]:
        # print(i)
        if i != "":
            sub=define_sub(i)
            sub_list.append(sub)
            gene_list.append(i)
        else:
            sub_list.append("")

    sub_list=list(filter(None,sub_list))
    # print(len(set(sub_list)))
    if len(set(sub_list)) == 1:
        pick+=1
        node.name=str(sub_list[0])+"_"+str(len(sub_list))+"_"+str(pick)
        collapse_gene_dic[node.name]=gene_list
        # print(node)
        return True
    else:
        # node.detach()
        return False
'''
相同科属种存在复杂拓扑结构,循环两次,collapse所有同一科属基因
'''
def collapsed_leaf2(node:Tree) -> bool:
    global pick2
    sub_list=[]
    gene_list=[]
    count=0
    for i in node2labels2[node]:
        # print(i)
        if i != "":
            sub=i.split("_")[0]
            num=int(i.split("_")[1])
            sub_list.append(sub)
            count+=num
            gene_list=gene_list+collapse_gene_dic[i]
        else:
            sub_list.append("")

    sub_list=list(filter(None,sub_list))
    # print(len(set(sub_list)))
    if len(set(sub_list)) == 1:
        pick2+=1
        node.name=str(sub_list[0])+"_"+str(count)+"_"+str(pick2)
        collapse_gene_dic2[node.name]=gene_list
        # print(node)
        return True
    else:
        # node.detach()
        return False


def main(tree_file:str,sup_value:int):
    gene_tree=Tree(tree_file)
    #V5增加长枝过滤
    
    '''V7
    引入重定根模块,获得尽可能多的clade HGT信息,但是运行时间增加多倍,未实装
    '''
    # gene_tree_list=re_root(gene_tree)
    
    outleaf=find_outlength(gene_tree)
    logging.info('{0} was delete due to extreme branch length'.format(" ".join(outleaf)))
    for i in gene_tree:
        if define_sub(i.name)=="OUT" or i.name in outleaf:
            i.delete()
    global node2labels1
    global collapse_gene_dic
    global pick
    node2labels1 = gene_tree.get_cached_content(store_attr="name")
    collapse_gene_dic={}
    pick=0

    # og=re.search('OG[0-9]{7}',tree_file).group()
    collapse_tree = Tree(gene_tree.write(is_leaf_fn=collapsed_leaf))
    #这里是关键的筛选步骤，1代表过滤掉在树结构中单独分布(1)的subgene


    for i in collapse_tree:
        num=int(i.name.split("_")[1])
        if num <= 1:
            i.delete()
    # print(collapse_tree)
    global node2labels2
    global collapse_gene_dic2
    global pick2
    node2labels2 = collapse_tree.get_cached_content(store_attr="name")
    collapse_gene_dic2={}
    pick2=0

    collapse_tree2 = Tree(collapse_tree.write(is_leaf_fn=collapsed_leaf2))
    # print(collapse_tree2)

    # find_clade(collapse_tree2)
    HGT,sus=find_clade(collapse_tree2,sup_value)
    if HGT:
        HGT_dic={}
        HGT_len_dic={}
        for i in HGT:
            if (collapse_tree2&i).up.support > sup_value:
                # print(og,True,collapse_gene_dic2[i],len(collapse_gene_dic2[i]),sep="\t")
                HGT_dic[i]=collapse_gene_dic2[i]
                HGT_len_dic[i]=len(collapse_gene_dic2[i])
        return HGT_dic,HGT_len_dic,True
    else:
        # print(og,False,"No HGT",0,sep="\t")
        return None,0,False
    if sus:
        for i in sus:
            if (collapse_tree2&i).up.support > sup_value:
                print(og,"sus",collapse_gene_dic2[i],len(collapse_gene_dic2[i]),sep="\t")    
if __name__ == "__main__":
    # tree_file=sys.argv[1]
    # tree_file=r"D:\study resource\HGT\5_find_clade\OG0000577_ali.mafft_rmdup_trimal_JTT.treefile"
    tree_file=r'D:\study resource\HGT\4_significant_gene\OG00003076_ali.mafft_rmdup_trimal_JTT.treefile'
    # sup_value=sys.argv[2]
    sup_value=50
    main(tree_file,sup_value)

# ###亚科列表
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
# OUT=["Plat","Acom"]

# ###亚科颜色列表
# color_template={"BB":"#b71515","OO":"#e97a0c","PP":"#ffde0a","PAC":"#034732","MAD":"#092e86","OUT":"#a9a29c"}