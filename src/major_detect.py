'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-11 20:55:49
LastEditors: Ne0tea
LastEditTime: 2025-02-10 21:49:11
'''
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
    V8 Python 3.7+
Author: Ne0tea
version: v8.0.0
Date: 2023-07-08 20:46:51
LastEditors: Ne0tea
LastEditTime: 2023-07-19 22:35:27
'''
import src.para_base as pb
import src.global_variable as glv
import sys
import time
from src.tree_top_parse.filter_length import find_outlength
from src.tree_top_parse.Inconsistant_toplogy_identify import decide_toplogy_from_process
from ete3 import Tree
import math

sys.setrecursionlimit(10000)
'''
根据gene_name确定科属 
INPUT:gene_name:str
OUTPUT:sub_family:str
'''
def get_subfamily_list(tree:Tree):
    root = tree
    if len(root.children) == 2:
        # 获取两个分组
        group1 = root.children[0]
        group2 = root.children[1]

        # 提取两个分组的叶子节点
        subfamily1 = [leaf.name for leaf in group1.iter_leaves()]
        subfamily2 = [leaf.name for leaf in group2.iter_leaves()]
    return [subfamily1, subfamily2]

def get_subfamily_list_from_phy(temp_phylo,OUT):
    if not OUT: #if there is no OUT in sub phylogeny, just return tree
        pass
    else:
        try:
            new_root = temp_phylo&OUT
            temp_phylo.set_outgroup(new_root)
            sister_root = new_root.get_sisters()[0]
            node=temp_phylo&OUT
            node.delete()
            temp_phylo=sister_root
        except:
            return get_subfamily_list(temp_phylo)
    return get_subfamily_list(temp_phylo)

def define_sub(gene):
    clade_subg_dic=glv.get('clade_subg_dic')
    sp=gene.split("_")[0]
    for i in clade_subg_dic:
        if sp in clade_subg_dic[i]:
            return i

'''
树结构,分布在相同拓扑结构中的sub_gene合并为同一枝
'''
def collapsed_leaf(node:Tree) -> bool:
    global pick
    sub_list=[]
    sp_list=[]
    gene_list=[]
    for i in node2labels1[node]:
        if i != "":
            sub=define_sub(i)
            sub_list.append(sub)
            sp_list.append(i.split("_")[0])
            gene_list.append(i)
        else:
            sub_list.append("")

    sub_list=list(filter(None,sub_list))
    sp_list=list(set(list(filter(None,sp_list))))
    if len(set(sub_list)) == 1:
        pick+=1
        ### Should be aligned to collapse format
        node.name=str(sub_list[0])+"_"+str('|'.join(sp_list))+"_"+str(pick)
        collapse_gene_dic[node.name]=gene_list
        return True
    else:
        return False

"""
将字符串每4个字符分割成一个列表
:param s: 输入字符串
:return: 由每4个字符组成的列表
"""
def split_string_every_4_chars(s):
    ### Should be aligned to collapse format
    cur_list=s.split("|")
    cur_list = [x for x in cur_list if x]
    # cur_list=[s[i:i+4] for i in range(0, len(s), 4)]
    return len(list(set(cur_list)))

def decide_out_name(clade_subg_dic):
    if 'OUT' in clade_subg_dic:
        return 'OUT'
    elif 'Out' in clade_subg_dic:
        return 'Out'
    elif 'Outgroup' in clade_subg_dic:
        return 'Outgroup'
    elif 'outgroup' in clade_subg_dic:
        return 'outgroup'
    else:
        return None

def read_config_phy(sub_phlo,OUT):
    if not OUT: #if there is no OUT in sub phylogeny, just return tree
        pass
    else:
        new_root = sub_phlo&OUT
        sub_phlo.set_outgroup(new_root)
        sister_root = new_root.get_sisters()[0]
        node=sub_phlo&OUT
        node.delete()
        sub_phlo=sister_root
    return sub_phlo

'''
根据基因树中的外类群对基因树重定根
INPUT:gene_tree:TREE
OUT:TREE:list 重定根后的树对象列表
'''
def re_root(gene_tree:Tree):
    outsp=[i for i in gene_tree.get_leaf_names() if define_sub(i) == "OUT" ]
    reatain_sp=gene_tree.get_leaf_names()
    reroot_tree_list=[]
    ascII_tree=gene_tree.write()
    for i in outsp:
        cur_tree=Tree(ascII_tree)
        cur_tree.set_outgroup(cur_tree&i)
        cur_tree.prune(reatain_sp)
        reroot_tree_list.append(cur_tree)

    if not outsp:
        cur_tree=Tree(ascII_tree)
        cur_tree.get_midpoint_outgroup()
        cur_tree.prune(reatain_sp)
        reroot_tree_list.append(cur_tree)
    return reroot_tree_list

'''
get leaf info of the sister node of current node
'''
def caluculate_sister_sp_componsent(node:Tree) -> dict:
    clade_dic={}
    cur_sister_nodes=node.get_sisters()
    for i in cur_sister_nodes:
        for leaf in i.get_leaf_names():
            if not leaf:
                continue
            clade=leaf.split("_")[0]
            number=leaf.split("_")[1]
            if clade in clade_dic:
                ### Should be aligned to collapse format
                clade_dic[clade]+='|'
                clade_dic[clade]+=number
            else:
                clade_dic[clade]=number
    return clade_dic

"""
判断 dict1 中是否存在任何一个键的值小于 dict2 中对应键的值
:param dict1: 第一个字典
:param dict2: 第二个字典
:return: 如果存在任意一个键满足 dict1[key] < dict2[key] 则返回 True，否则返回 False
"""
def any_key_less(dict1, dict2):
    for key, value in dict2.items():
        if key in dict1 :
            dict1_vlaue=split_string_every_4_chars(dict1[key])
            if value > dict1_vlaue:
                return True
        else:
            return True
    return False

def all_key_less(dict1, dict2):
    for key, value in dict2.items():
        if key in dict1 :
            dict1_vlaue=split_string_every_4_chars(dict1[key])
            if value <= dict1_vlaue:
                return False
    return True

"""
判断 dict1 中是否存在任何一个键的值小于 dict2 中对应键的值
:param dict1: 第一个字典
:param dict2: 第二个字典
:return: 如果存在任意一个键满足 dict1[key] < dict2[key] 则返回 True，否则返回 False
"""
def less_half_key_overlap(dict1, dict2, template_dic):
    high_keys1 = {k: v for k, v in dict1.items() if k in template_dic and split_string_every_4_chars(v) > template_dic[k]}
    high_keys2 = {k: v for k, v in dict2.items() if k in template_dic and split_string_every_4_chars(v) > template_dic[k]}

    template_half_count = math.ceil((len(template_dic)-1) / 2)

    high_keys1_set = set(high_keys1.keys())
    high_keys2_set = set(high_keys2.keys())

    if high_keys2_set <= high_keys1_set and len(high_keys2_set) >= template_half_count:
        return False
    else:
        return True

def remove_duplicate_matrices(matrix_dict):
    """去除字典中重复的 NumPy 矩阵值，返回新的字典。"""
    seen = set()
    new_dict = {}
    for key, matrix in matrix_dict.items():
        matrix_tuple = tuple(map(tuple, matrix))
        if matrix_tuple not in seen:
            new_dict[key] = matrix
            seen.add(matrix_tuple)
    return new_dict

def get_ample_step(last_ample_clade,dict3,temp_dic):
    out_key=last_ample_clade.copy()
    ample_key2=[]
    for key in dict3:
        vlaue=split_string_every_4_chars(dict3[key])
        if key in temp_dic and vlaue >= temp_dic[key]:
            ample_key2.append(key)
    out_key.extend(ample_key2)

    return out_key

def exist_sister_node(node,node_list, subfamily_list):
    cur_node_list=node_list.copy()
    if node not in cur_node_list:
        return False
    cur_node_list.remove(node)
    node_name = node.name
    if node_name.split('_')[0] in subfamily_list[0]:
        cur_subf = 'subf1'
    elif node_name.split('_')[0] in subfamily_list[1]:
        cur_subf = 'subf2'
    else:
        return False
    for i in cur_node_list:
        now_name = i.name
        if now_name.split('_')[0] in subfamily_list[0]:
            now_subf = 'subf1'
        elif now_name.split('_')[0] in subfamily_list[1]:
            now_subf = 'subf2'
        else:
            continue
        if node.up == i.up and now_subf == cur_subf:
            return True
    return False

"""
根据叶子节点获取进化过程

参数:
internal_node: Tree类型, 表示树中的一个节点
Clade_freq_dic_template: dict类型, 模板克隆频率字典
Clade_template_dic: dict类型, 模板克隆信息字典

返回值:
cur_ample_process: 列表, 当前节点及其祖先节点的克隆频率达到template要求时的进化过程,根据iter顺序添加
cur_node_clade_freq: 字典, 当前节点的克隆频率
iter_node_list: 列表, 进化过程中遍历的节点列表
og_top_node: SG group 的顶层节点
"""
def get_process_from_leaf(internal_node:Tree,Clade_freq_dic_template,Clade_template_dic):
    cur_clade=internal_node.name.split('_')[0]
    cur_number=internal_node.name.split('_')[1]##!!!notice this cur number is DbraDbraBamp (str)

    cur_node_clade_freq={x: '' for x in Clade_freq_dic_template}
    cur_node_clade_freq[cur_clade]+=cur_number

    cur_ample_process=[[x for x in cur_node_clade_freq if split_string_every_4_chars(cur_node_clade_freq[x]) >= Clade_template_dic[x] ]]
    cur_clade_process=[[cur_clade]] # record the clade within the iteral process
    #self -> sister -> up/sister -> up/sister ...
    cur_node=internal_node
    iter_node_list=[cur_node] # record the iter node within the genetree
    CSC=caluculate_sister_sp_componsent(cur_node)
    og_top_node=internal_node
    while any_key_less(cur_node_clade_freq,Clade_template_dic) and any_key_less(CSC,Clade_template_dic) and \
        less_half_key_overlap(cur_node_clade_freq, CSC, Clade_template_dic) and cur_node.up:
        # if 'PAN' in internal_node.name and 'Msin−B' in internal_node.name and 'Sspo−A' in internal_node.name and \
        #     'Sbic' in internal_node.name:
        #         print(cur_node_clade_freq)
        cur_ample_step=get_ample_step(cur_ample_process[-1],CSC,Clade_template_dic)
        if cur_ample_step == cur_ample_process[-1]:
            last_CSC=CSC.copy()
        else:
            last_CSC={}
        cur_ample_process.append(cur_ample_step)
        cur_clade_process.append([x for x in CSC if CSC[x]])
        for i in CSC:
            if i in cur_node_clade_freq:
                if cur_node_clade_freq[i]:
                    cur_node_clade_freq[i]+='|'
                cur_node_clade_freq[i]+=CSC[i]
            else:
                cur_node_clade_freq[i]=CSC[i]
        # if 'BB' in internal_node.name and 'Phedu' in internal_node.name and 'DlatC' in internal_node.name and 'Bamp' in internal_node.name:
        #     print(cur_node_clade_freq)
        cur_node=cur_node.up
        og_top_node=og_top_node.up
        iter_node_list.append(cur_node)
        CSC=caluculate_sister_sp_componsent(cur_node)
        # if 'BB' in internal_node.name and 'Phedu' in internal_node.name and 'DlatC' in internal_node.name and 'Bamp' in internal_node.name:
        #     print(CSC)
        if last_CSC:
            for i in last_CSC:
                if i in CSC:
                    if CSC[i]:
                        CSC[i]+='|'
                    CSC[i]+=last_CSC[i]
                else:
                    CSC[i]=last_CSC[i]

    return cur_ample_process,cur_node_clade_freq,iter_node_list,og_top_node,cur_clade_process

def get_contary_leaf(tree_node:Tree,subfamily):
    cur_node_agf_gene_list=[]
    for leaf in tree_node.get_leaf_names():
        if leaf.split("_")[0] not in subfamily:
            cur_node_agf_gene_list.append(leaf)
    return cur_node_agf_gene_list

"""
iteral SG tree
:param cur_node: 当前的节点对象
:param Clade_template_dic: 模板字典，用于比对
:param support_value: support value of the branch
:return: 布尔值，donor clade，donor sp, receive clade, receive sp
"""
def judge_sg_leaf(cur_node:Tree,Clade_template_dic,subfamily_list,support_value:int):
    cur_template_dic={x:Clade_template_dic[x] for x in subfamily_list[0]+subfamily_list[1]}
    subf1, subf2 = subfamily_list
    cur_clade = cur_node.name.split('_')[0]
    cur_subf1 = 'subf1' if cur_clade in subf1 else 'subf2'
    sister_clade_count:dict=caluculate_sister_sp_componsent(cur_node)
    # if cur_node.name.split('_')[0] == 'CHL' and cur_node.name.split('_')[1] == 'Cson-B':
    #     print(cur_node.up)
    if not cur_node.up:
        return False, None, None, None, None, None, None
    if cur_node.up.support < support_value:
        return False, None, None, None, None, None, None
    node_for_clade = cur_node
    node_for_sp = cur_node

    if all_key_less(sister_clade_count,cur_template_dic):
        # sister node didn't contain sample clade
        while node_for_clade.up and all_key_less(sister_clade_count,cur_template_dic):
            node_for_clade = node_for_clade.up
            sister_clade_count=caluculate_sister_sp_componsent(node_for_clade)
        if all_key_less(sister_clade_count,cur_template_dic):
            return False, None, None, None, None, None, None
        ample_clade = get_ample_step([],sister_clade_count,cur_template_dic)
        if len(ample_clade) > 1 or ample_clade[0] == cur_clade:
            return False, None, None, None, None, None, None
        else:
            donor_clade = ample_clade[0]
    else:
        # sister node is regared as Donor clade
        ample_clade = get_ample_step([],sister_clade_count,cur_template_dic)
        if len(ample_clade) > 1 or ample_clade[0] == cur_clade:
            return False, None, None, None, None, None, None
        else:
            donor_clade = ample_clade[0]

    # if Donor subf same as cur subf False
    donor_subf = 'subf1' if donor_clade in subf1 else 'subf2'
    if donor_subf == cur_subf1:
        return False, None, None, None, None, None, None

    # node up until the node contain the same clade
    sister_clade_list=[]
    for sister_node in node_for_sp.get_sisters():
        sister_clade_list.extend(sister_node.get_leaf_names())
    while node_for_sp.up and donor_clade not in [x.split('_')[0] for x in sister_clade_list]:
        node_for_sp = node_for_sp.up
        sister_clade_list=[]
        for sister_node in node_for_sp.get_sisters():
            sister_clade_list.extend(sister_node.get_leaf_names())

    # Check if the sister clade is the same as the donor clade
    # Check node sister consistant to the donor clade
    if [donor_clade] != [x.split('_')[0] for x in sister_clade_list]:
        return False, None, None, None, None, None, None
    else:
        node_for_sp = node_for_sp.up
        sister_sp = []
        for i in node_for_sp.get_sisters():
            sister_sp.extend(i.get_leaf_names())
        if donor_clade in [x.split('_')[0] for x in sister_sp]:
            donor_sp = '|'.join([x.split('_')[1] for x in sister_sp if x.split('_')[0] == donor_clade])
            receive_clade = cur_clade
            receive_sp = cur_node.name.split('_')[1]
            receive_gene = collapse_gene_dic[cur_node.name]
            donor_gene = []
            for donor_collapse in [x for x in sister_sp if x.split('_')[0] == donor_clade]:
                donor_gene.extend(collapse_gene_dic[donor_collapse])
            donor_gene = '|'.join(donor_gene)
            return True, donor_clade, donor_sp, receive_clade, receive_sp, donor_gene, receive_gene
        else:
            return False, None, None, None, None, None, None

'''
根据collapsed基因树,过滤长枝之后,鉴定AGF产生的与物种树拓扑inconstant的高bootstrap值的Gene Set:
INPUT:collapsed_gene_tree:Tree
      support_value:int
OUTPUT:return HGT_gene_list number_of_HGT_gene
'''
def find_inconstant_node(collapsed_gene_tree:Tree, rootnested, out_name, uncollapsed_tree:Tree,sup_value=50, test_mode:bool=False):
    sub_phlo=glv.get('sub_phlo')
    subfamily_list=get_subfamily_list_from_phy(Tree(sub_phlo),out_name)
    clade_subg_dic=glv.get("clade_subg_dic")

    # Clade_template_dic set the createria of clade number in subtree
    Clade_template_dic={x:int(math.sqrt(len(clade_subg_dic[x]))) for x in clade_subg_dic if x != out_name }
    Clade_template_dic[out_name]=1
    # print(Clade_template_dic)
    '''
    {'PANCHL': ['PAN', 'CHL'], 'POOORYBAM': ['POO', 'ORY', 'BAM']}
    '''
    Clade_freq_dic_template={x:'' for x in clade_subg_dic}
    AGF_gene_set=[]
    same_clade_in_AGF_set = []
    All_GF_gene_node=[]
    AGF_gene_node=[]
    SG_top_node=[]
    if rootnested:
        NAGF_gene_set=[]
        NSG_top_node=[]
    for internal_node in collapsed_gene_tree:

        # Skip node contain less species than template
        cur_unit=internal_node.name.split('_')
        if cur_unit[0]==out_name or split_string_every_4_chars(cur_unit[1]) < Clade_template_dic[cur_unit[0]]:
            continue

        #plus current clade freq into dic
        cur_ample_process,cur_node_clade_freq,iter_node_list,top_node,cur_clade_process=get_process_from_leaf(internal_node,Clade_freq_dic_template,Clade_template_dic)

        for clade in cur_node_clade_freq:
            vlaue=split_string_every_4_chars(cur_node_clade_freq[clade])
            cur_node_clade_freq[clade]=vlaue
        # print(cur_ample_process)
        # if 'PAN' in internal_node.name and 'Cmac' in internal_node.name and 'Pgig' in internal_node.name and \
        #     'Zmay' in internal_node.name:
        #     toplogy,cur_clade_name,clade_index=decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name)
        #     print(cur_ample_process,clade_index)
        #     if toplogy != 'Outlier_clade' and toplogy != 'Inconstant_clade' and toplogy != 'Norooted':
        #         print(cur_ample_process,toplogy,top_node)
        #         print(subfamily_list)
        #         print(decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name))
        #         print(iter_node_list[clade_index])
        if len(cur_ample_process) >= 2:
            if out_name in cur_ample_process[-2]:
                if rootnested:
                    for susp in top_node:
                        if susp.name.split("_")[0] == out_name:
                            top_node.set_outgroup(susp)
                            for leaf in top_node:
                                cur_ample_process,cur_node_clade_freq,iter_node_list,_,cur_clade_process=get_process_from_leaf(leaf,Clade_freq_dic_template,Clade_template_dic)
                                toplogy,cur_clade_name,clade_index=decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name)
                                if toplogy != 'Outlier_clade' and toplogy != 'Inconstant_clade':
                                    cur_agf_gene_list=get_contary_leaf(iter_node_list[clade_index],cur_clade_name)
                                    NAGF_gene_set.extend(cur_agf_gene_list)
                            NSG_top_node.append(top_node)
                else:
                    continue
            else:
                toplogy,cur_clade_name,clade_index=decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name,cur_clade_process)
        else:
            toplogy,cur_clade_name,clade_index=decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name,cur_clade_process)
        subf_sp_list=subfamily_list[0] if cur_clade_name == 'subf1' else subfamily_list[1]

        ###toplogy can be AGF, TGF, Outlier_clade, Inconstant_clade, or Norooted
        if toplogy != 'Outlier_clade' and toplogy != 'Inconstant_clade' and toplogy != 'Norooted':
            # cur_agf_contary_gene_list=[]
            if toplogy == 'AGF':
                AGF_gene_node.append(iter_node_list[clade_index])
                cur_agf_gene_list=get_contary_leaf(iter_node_list[clade_index],subf_sp_list)
                # cur_agf_gene_list, cur_agf_contary_gene_list = get_contary_leaf_in_AGF(iter_node_list[clade_index],subf_sp_list)
            elif toplogy == 'TGF':
                cur_agf_gene_list=get_contary_leaf(iter_node_list[clade_index],subf_sp_list)

            All_GF_gene_node.append(iter_node_list[clade_index])
            AGF_gene_set.append(cur_agf_gene_list)
            SG_top_node.append(top_node)

    rGF_dic={}
    for sg_tree in SG_top_node:
        for nnnn in sg_tree:
            if nnnn.name.split("_")[0] == out_name:
                sg_tree.set_outgroup(nnnn)
                break
        all_leaf = []
        for x in sg_tree.get_leaf_names():
            all_leaf+=collapse_gene_dic[x]
        uncollapsed_tree_node = uncollapsed_tree.get_common_ancestor(all_leaf)
        sg_tree_newick=uncollapsed_tree_node.write()
        sg_tree_node_for_iter = sg_tree
        for sg_leaf in sg_tree_node_for_iter:
            cur_unit=sg_leaf.name.split('_')
            if cur_unit[0]==out_name or split_string_every_4_chars(cur_unit[1]) >= Clade_template_dic[cur_unit[0]]:
                continue
            status, sg_donor_clade, sg_donor_sp, sg_receive_clade, sg_receive_sp, sg_donor_gene, sg_receive_gene = judge_sg_leaf(sg_leaf,Clade_template_dic,subfamily_list,sup_value)

            if status:
                rGF_dic[sg_leaf.name] = (sg_donor_clade, sg_donor_sp, sg_receive_clade, sg_receive_sp, sg_donor_gene, sg_receive_gene, sg_tree_newick)

    if AGF_gene_set:
        ###format [['PP_Aspl|Telo_7', 'BB_DlatC_9', 'BB_Phedu_11', 'BB_Bamp|DlatC_13', 'BB_Olat_15'],
        ###         ['PAC_Zmay|MsinB|Asem|Sbic|Doli|EcruB|EcolD|Sita|EcruC|EcruA|EoryB|Ehap|EcolF|DexiB|MsinA|Phal|EcolE|PvirN|EoryA|Svir_189'], 
        ###         ['MAD_CsonB|Otho|Ecur|EtefB|Enin|CsonA|Zjap|EtefA_191'], [], []]
        AGF_dic={}
        AGF_len_dic={}
        AGF_SG_node={}
        ### SUPPORT Value was set as creatira of gene flow node
        for idx,collapse_node in enumerate(All_GF_gene_node):
            top_node_all_sp=[]
            if not collapse_node:
                continue
            if exist_sister_node(collapse_node, AGF_gene_node, subfamily_list):
                if collapse_node.up.up.support > sup_value:
                    for collapse_leaf in AGF_gene_set[idx]:
                        AGF_dic[collapse_leaf]=collapse_gene_dic[collapse_leaf]
                        AGF_len_dic[collapse_leaf]=len(collapse_gene_dic[collapse_leaf])
                        for subsetleaf in SG_top_node[idx].get_leaf_names():
                            top_node_all_sp.extend(collapse_gene_dic[subsetleaf])
                        uncollapsed_node = uncollapsed_tree.get_common_ancestor(top_node_all_sp)
                        AGF_SG_node[collapse_leaf]=uncollapsed_node
            else:
                if collapse_node.up.support > sup_value:
                    for collapse_leaf in AGF_gene_set[idx]:
                        if (collapsed_gene_tree&collapse_leaf).up.support > sup_value:
                            AGF_dic[collapse_leaf]=collapse_gene_dic[collapse_leaf]
                            AGF_len_dic[collapse_leaf]=len(collapse_gene_dic[collapse_leaf])
                        for subsetleaf in SG_top_node[idx].get_leaf_names():
                            top_node_all_sp.extend(collapse_gene_dic[subsetleaf])
                        uncollapsed_node = uncollapsed_tree.get_common_ancestor(top_node_all_sp)
                        AGF_SG_node[collapse_leaf]=uncollapsed_node

        if rootnested:
            return AGF_dic, AGF_len_dic, AGF_SG_node, list(set(NAGF_gene_set)), list(set(NSG_top_node)), rGF_dic
        else:
            return AGF_dic, AGF_len_dic, AGF_SG_node, rGF_dic
    else:
        return None, None, None, None

def main(tree_file:str,sup_value:int, rootnested,target_name=None,test_mode:bool=False):
    clade_subg_dic=glv.get("clade_subg_dic")
    out_name=decide_out_name(clade_subg_dic)
    gene_tree=Tree(tree_file)

    ##expand node support value if all node support value is less than 1
    all_support_less_than_1 = all(node.support <= 1 for node in gene_tree.traverse() if not node.is_leaf())
    if all_support_less_than_1:
        for node in gene_tree.traverse():
            if not node.is_leaf():
                node.support *= 100

    #leaf length filteration
    outleaf=find_outlength(gene_tree,clade_subg_dic[out_name])
    for i in gene_tree:
        if i.name in outleaf:
            i.delete()
    global node2labels1
    global collapse_gene_dic
    global pick
    node2labels1 = gene_tree.get_cached_content(store_attr="name")
    collapse_gene_dic={}
    pick=0

    AGF_dic={}
    AGF_len_dic={}
    AGF_top_node={}
    rGF_dic={}
    back_up_tree=Tree(gene_tree.write())
    collapse_tree = Tree(gene_tree.write(is_leaf_fn=collapsed_leaf))
    if rootnested:
        c_AGF_dic,c_AGF_len_dic,c_AGF_top_node,NAGF_gene_list,Ngenenode_set, c_rgf_dic=find_inconstant_node(collapse_tree,rootnested,out_name,back_up_tree,sup_value,test_mode)
        if c_AGF_dic:
            AGF_dic,AGF_len_dic,AGF_top_node, rGF_dic=c_AGF_dic,c_AGF_len_dic,c_AGF_top_node, c_rgf_dic
        re_rooted_tree=re_root(back_up_tree)
        for rtree in re_rooted_tree:
            # rtree_backup=Tree(rtree.write())
            pick=0
            collapse_gene_dic={}
            node2labels1 = rtree.get_cached_content(store_attr="name")
            cur_collapse_tree = Tree(rtree.write(is_leaf_fn=collapsed_leaf))
            c_AGF_dic,c_AGF_len_dic,c_AGF_top_node,NAGF_gene_list,Ngenenode_set, c_rgf_dic=find_inconstant_node(cur_collapse_tree,rootnested,out_name,back_up_tree,sup_value,test_mode)
            if c_AGF_dic:
                AGF_dic.update(c_AGF_dic)
                AGF_len_dic.update(c_AGF_len_dic)
                AGF_top_node.update(c_AGF_top_node)
            if c_rgf_dic:
                rGF_dic.update(c_rgf_dic)
    else:
        c_AGF_dic,c_AGF_len_dic,c_AGF_top_node, c_rgf_dic=find_inconstant_node(collapse_tree,rootnested,out_name,back_up_tree,sup_value,test_mode)
        if c_AGF_dic:
            AGF_dic,AGF_len_dic,AGF_top_node, rGF_dic=c_AGF_dic,c_AGF_len_dic,c_AGF_top_node, c_rgf_dic
        re_rooted_tree=re_root(back_up_tree)
        # For test
        # back_up_tree.set_outgroup(back_up_tree&'Plat_prot_Plat_Pl01g19720_mRNA1')
        # re_rooted_tree=[back_up_tree]
        # print(len(re_rooted_tree))
        for rtree in re_rooted_tree:
            # rtree_backup=Tree(rtree.write())
            pick=0
            collapse_gene_dic={}
            node2labels1 = rtree.get_cached_content(store_attr="name")
            cur_collapse_tree = Tree(rtree.write(is_leaf_fn=collapsed_leaf))

            c_AGF_dic,c_AGF_len_dic,c_AGF_top_node, c_rgf_dic=find_inconstant_node(cur_collapse_tree,rootnested,out_name,back_up_tree,sup_value,test_mode)
            if c_AGF_dic:
                AGF_dic.update(c_AGF_dic)
                AGF_len_dic.update(c_AGF_len_dic)
                AGF_top_node.update(c_AGF_top_node)
            if c_rgf_dic:
                rGF_dic.update(c_rgf_dic)

    unique_value_dict = {}
    seen_values = []
    for key, value in AGF_dic.items():
        if value not in seen_values:
            seen_values.append(value)
            unique_value_dict[key] = value
    AGF_len_dic={x: y for x, y in AGF_len_dic.items() if x in unique_value_dict}
    AGF_top_node={x: y for x, y in AGF_top_node.items() if x in unique_value_dict}

    unique_c_rgf_dic = {}
    seen_values = []
    for key, value in rGF_dic.items():
        if value[4] not in seen_values:
            seen_values.append(value[4])
            unique_c_rgf_dic[key] = value

    if unique_value_dict:
        if rootnested:
            return unique_value_dict,AGF_len_dic,AGF_top_node,True,Ngenenode_set, unique_c_rgf_dic
        else:
            return unique_value_dict,AGF_len_dic,AGF_top_node,True, unique_c_rgf_dic
    else:
        if rootnested:
            return None, None, None, False, Ngenenode_set, unique_c_rgf_dic
        else:
            return None, None, None, False, unique_c_rgf_dic

if __name__ == "__main__":
    tree_file=r'E:\Bio_analysis\HGT_newcollection\1_gene_tree_set\OG0000328.fa.rd.tml.JTT.contree'
    # tree_file=r'E:\Bio_analysis\HGT_newcollection\HGTfinder_script\AGFD_V2_0910\Example_file\Testfile3.treefile'
    sup_value=70
    config=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    # config=r'E:\Bio_analysis\HGT_newcollection\HGTfinder_script\AGFD_V2_0910\Example_file\test_config.txt'
    pb.main(config)
    clade_subg_dic=glv.get('clade_subg_dic')

    start_time=time.time()
    rootnested=False
    AGF_gene_dic, AGF_gene_number_dic, AGF_top_node_dic, AGF_status, rGF_dic=main(tree_file,sup_value,rootnested,'Pgi03A00007030',False)
    end_time=time.time()
    print(end_time-start_time)

    template_color={"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1",'Target':'lightgrey'}
    # for i in AGF_gene_dic:
    #     AGF_gene_list = AGF_gene_dic[i]
    #     # if  AGF_gene_number_dic[i] != 12: continue
    #     SG_node = AGF_top_node_dic[i]
    #     print(len(AGF_top_node_dic[i].get_leaf_names()))
    #     print(AGF_gene_list)
    #     donor_clade='BAM'
        # modified_branch_lenth_test.run_test(SG_node,AGF_gene_list,template_color,donor_clade,clade_subg_dic)
    for i in rGF_dic:
        donor_clade, donor_sp, receive_clade, receive_sp, donor_gene, receive_gene, sg_newick= rGF_dic[i]
        print(donor_clade, donor_sp, receive_clade, receive_sp, receive_gene)