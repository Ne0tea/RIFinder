'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-12-31 23:53:19
LastEditors: Ne0tea
LastEditTime: 2025-01-09 17:51:32
'''
import sys
import re
import numpy as np
from ete3 import Tree

genetree_file=sys.argv[1]
# genetree_file=r'E:\Bio_analysis\HGT_newcollection\GF_4295_sg_prot.contree'
gene_tree=Tree(genetree_file)
leaf_gene_tree_list = gene_tree.get_leaf_names()
out_gene = [x for x in leaf_gene_tree_list if x.split('_')[0] in ['Plat','Acom']]
if out_gene:
    gene_tree.set_outgroup(out_gene[0])

og_id=re.search('(GF_[0-9]*)_',genetree_file).group(1)

AGF_result_file=sys.argv[2]
# AGF_result_file=r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out.consistance.significant.all.out'
# gene1=sys.argv[2]
with open(AGF_result_file,'r') as f:
    for line in f:
        line=line.strip()
        c_id=line.split('\t')[0]
        if og_id == c_id:
            cur_receive_clade=line.split('\t')[5]
            cur_receive_sp = line.split('\t')[8]
            cur_receive_gene = line.split('\t')[9].split('|')
            cue_donor_sp = line.split('\t')[6]
cur_receive_gene = [re.search(r'_prot_(.*)', x).group(1) for x in cur_receive_gene]
cur_receive_gene = [x for x in cur_receive_gene if x in leaf_gene_tree_list]
if not cur_receive_gene:
    print("No available gene! exit")
    exit()
#outfile=gene1+"_hypothsis"
MCA_receive_node = gene_tree.get_common_ancestor(list(set(cur_receive_gene)))
MCA_receive_gene_list = MCA_receive_node.get_leaf_names()

BB=["Aluo-C", "Aluo-D", "Bamp-A", "Bamp-B", "Bamp-C", "Dbra1-A", "Dbra1-B", "Dbra1-C", "Dbra2-A", "Dbra2-B", "Dbra2-C", "Dlat-A", "Dlat-B", "Dlat-C", "Dsin-A", "Dsin-B", "Dsin-C", "Gang-B", "Gang-C", "Hcal-C", "Hcal-D", "Mbac-A", "Mbac-B", "Mbac-C", "Olat", "Otgla-B", "Otgla-C", "Pedu-C", "Pedu-D", "Rgui", "Rrac-B", "Rrac-C"]
OO=["Lper", "Oalt-C", "Oalt-D", "Obar", "Obra", "Ogla", "Oglu", "Omer", "Opun", "Oruf", "Osat", "Zlat"]
PP=["Aatl", "Aeri", "Ains-C", "Ains-D", "Alon", "Amyo", "Aspl", "Astr", "Atau", "Aumb", "Bdis", "Bsta", "Dvil", "Hmar", "Hvul-C", "Hvul-S", "Lmul", "Pten", "Scer", "Taes-A", "Taes-B", "Taes-D", "Tdic-A", "Tdic-B", "Telo", "Tmon", "Ttur-A", "Ttur-B", "Tura"]
PAN=["Asem", "Came", "Cmac-A", "Cmac-B", "Cpur-A", "Cpur-B", "Dexi-A", "Dexi-B", "Dsan-C", "Dsan-D", "Dsan-E", "Ecol-D", "Ecol-E", "Ecol-F", "Ecru-A", "Ecru-B", "Ecru-C", "Ehap", "Eoph", "Eory-A", "Eory-B", "Msin-A", "Msin-B", "Pgig-A", "Pgig-B", "Phal", "Ppur-A", "Ppur-B", "Pvag", "Pvir-K", "Pvir-N", "Sbic", "Sita", "Sspo-A", "Sspo-B", "Sspo-C", "Sspo-D", "Svir", "Zmay"]
CHL=["Ecur", "Enin", "Cson-A", "Cson-B", "Eind", "Etef-A", "Etef-B", "Otho", "Zjap"]
OUT=["Plat","Acom"]

def define_sub(gene):
    sp=gene.split("_")[0]
    if sp in BB:
        return "BAM"
    elif sp in OO:
        return "ORY"
    elif sp in PP:
        return "POO"
    elif sp in PAN:
        return "PAN"
    elif sp in CHL:
        return "CHL"
    elif sp in OUT:
        return "OUT"
    else:
        return "empty"

GF_node_list = []
opp_dic={}
fam=cur_receive_clade
#store leaves belong to the same subfamily with transfered genes
GF_status = 'ILS'
for i in gene_tree.iter_leaf_names():
    c_sp=i.split('_')[0]
    c_fam=define_sub(c_sp)
    if c_fam==fam and i not in MCA_receive_gene_list:
        GF_status = 'GF'
        GF_node_list.append(i)
    elif c_fam!=fam:
        if c_fam in opp_dic:
            opp_dic[c_fam].append(i)
        else:
            opp_dic[c_fam]=[i]
if GF_status == 'GF':
    GF_node_dis = [(x, np.mean([(gene_tree&x).get_distance(y) for y in cur_receive_gene])) for x in GF_node_list]
else:
    if cur_receive_clade == 'BAM':
        related_clade = "POO"
    elif cur_receive_clade == 'ORY':
        related_clade = "BAM"
    elif cur_receive_clade == 'POO':
        related_clade = "BAM"
    elif cur_receive_clade == 'PAN':
        related_clade = "CHL"
    elif cur_receive_clade == 'CHL':
        related_clade = "PAN"
    related_gene_list = opp_dic[related_clade]
    GF_node_dis = [(x, np.mean([(gene_tree&x).get_distance(y) for y in cur_receive_gene])) for x in related_gene_list]
insert_gene_name=min(GF_node_dis,key=lambda x:x[1])[0]

#Move the transfer gene to hyposized position
MCA_receive_node_backup = MCA_receive_node.copy()
center=gene_tree&insert_gene_name
removed_node = MCA_receive_node.detach()

R_node = center.add_sister(name='Receive_node')
R_node.add_child(MCA_receive_node_backup)


gene_tree.prune(list(set(cur_receive_gene)), preserve_branch_length=True)
gene_tree.write(format=1, outfile=og_id+"_hypothsis")