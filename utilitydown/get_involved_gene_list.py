'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-29 16:56:20
LastEditors: Ne0tea
LastEditTime: 2024-11-07 11:02:02
'''
import re
from ete3 import Tree
import time
import argparse
template_color={"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1"}

def timer(func):
    def func_wrapper(*args, **kwargs):
        from time import time
        time_start = time()
        result = func(*args, **kwargs)
        time_end = time()
        time_spend = time_end - time_start
        print('%s cost time: %.3f s' % (func.__name__, time_spend))
        return result
    return func_wrapper
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

def make_config(config_file):
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
        if x != decide_out_name(clade_subg_dic):
            sp_ord.extend(clade_subg_dic[x])
    return clade_subg_dic,subg_clade_dic,sp_ord

@timer
def main(ra_out_file, rgf_out_file, ra_nwk_file, config_file):
    clade_subg_dic,subg_clade_dic,sp_ord = make_config(config_file)

    ###Build GF sg tree dic
    nwk_dic={}
    with open(ra_nwk_file,'r') as nwkf:
        for i in nwkf:
            line=i.strip().split()
            nwk_dic[line[0]]=line[3]


    out_file = open(ra_out_file+'.gene.list','w')
    with open(ra_out_file,'r') as ra_out:
        for ral in ra_out:
            line=ral.strip().split('\t')
            '''
            out file format example:
            donor_clade, receive_clade, donor_sp, receive_sp
            GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
            '''
            donor_clade, receive_clade, donor_sp, donor_gene, receive_sp, receive_gene= line[4], line[5],line[6].split('|'), line[7].split('|'), line[8].split('|'), line[9].split('|')
            tree = Tree(nwk_dic[line[0]])
            native_gene_dict = {}
            for i in tree.get_leaf_names():
                # gene_clade=subg_clade_dic[i.split('_')[0]]
                gene_sp = i.split('_')[0]
                if i in receive_gene or i in donor_gene : continue
                if gene_sp in receive_sp:
                    i = re.search(r'_prot_(.*)', i).group(1)
                    if gene_sp in native_gene_dict:
                        native_gene_dict[gene_sp].append(i)
                    else:
                        native_gene_dict[gene_sp]=[i]

            for i in receive_gene:
                if 'prot' in i:
                    i = re.search(r'_prot_(.*)', i).group(1)
                cur_line = '\t'.join([line[0], 'receive', receive_clade, i])
                out_file.write(cur_line+'\n')
            for sp in native_gene_dict:
                native_gene_list = native_gene_dict[sp]
                for ng in native_gene_list:
                    cur_line = '\t'.join([line[0], 'native', receive_clade, ng])
                    out_file.write(cur_line+'\n')
            for i in donor_gene:
                if 'prot' in i:
                    i = re.search(r'_prot_(.*)', i).group(1)
                cur_line = '\t'.join([line[0], 'donor', receive_clade, i])
                out_file.write(cur_line+'\n')

    with open(rgf_out_file,'r') as ra_out:
        for ral in ra_out:
            line=ral.strip().split('\t')
            '''
            out file format example:
            donor_clade, receive_clade, donor_sp, receive_sp
            GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
            '''
            trans_type, donor_clade, receive_clade, donor_sp, donor_gene, receive_sp, receive_gene= line[3], line[4], line[5],line[6].split('|'), line[7].split('|'), line[8].split('|'), line[9].split('|')
            if trans_type != 'RGF': continue
            for i in receive_gene:
                if 'prot' in i:
                    i = re.search(r'_prot_(.*)', i).group(1)
                cur_line = '\t'.join([line[0], 'receive', receive_clade, i])
                out_file.write(cur_line+'\n')
            for sp in native_gene_dict:
                native_gene_list = native_gene_dict[sp]
                for ng in native_gene_list:
                    cur_line = '\t'.join([line[0], 'native', receive_clade, ng])
                    out_file.write(cur_line+'\n')
            for i in donor_gene:
                if 'prot' in i:
                    i = re.search(r'_prot_(.*)', i).group(1)
                cur_line = '\t'.join([line[0], 'donor', receive_clade, i])
                out_file.write(cur_line+'\n')
    out_file.close()

if __name__ == "__main__":
    begin = time.perf_counter()
    # parser = argparse.ArgumentParser(description="Stat donor/receiver from the AGFD result.")
    # parser.add_argument("-i","--outfile", type=str, help="Outfile from AGFD.")
    # parser.add_argument("-r","--rawout", type=str, help="raw Outfile from AGFD.")
    # parser.add_argument("-t","--outtree", type=str, help="Outfile of SG tree from AGFD.")
    # parser.add_argument("-c","--config", type=str, help="Path to the config file.")
    # args = parser.parse_args()
    # main(args.outfile, args.rawout, args.outtree, args.config)
    '''
    out file format example:
    donor_clade, receive_clade, donor_sp, receive_sp
    OG0000020	Sbic_prot_Sbic_Sobic_004G235600_1_v3_2|Eoph_prot_Eoph_ctg632_38_m	2   ORY	PAN	Ogla|Oruf|Osat	Eoph|Sbic
    '''
    ra_out_file=r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out.consistance.significant'
    rgf_out_file=r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out'
    ra_nwk_file = r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out.OG.node.nwk'
    config_file=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'

    main(ra_out_file, rgf_out_file, ra_nwk_file, config_file)

    end = time.perf_counter()
    h=int((end - begin)/3600)
    m=int((end - begin)/60 - h*60)
    s=int(end - begin - h*3600 - m*60)
    print("All tasks used time: %sh %sm %ss" % (h,m,s))