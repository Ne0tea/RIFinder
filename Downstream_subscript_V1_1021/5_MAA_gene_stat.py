'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-09-19 23:40:14
LastEditors: Ne0tea
LastEditTime: 2023-10-15 19:04:39
'''

import parse_gff3
from ete3 import Tree
import pandas as pd
import numpy as np
import difflib
import os
import sys
import re
import multiprocessing
from datetime import datetime

import traceback

def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result
    pass

def get_gff_gene_base(gff_list:list,gene_size=str):
    t_intro={}
    g_size={}
    for i in gff_list:
        print(os.path.basename(i),'done')
        sp=os.path.basename(i.strip()).split('_')[0]
        try:
            t_intro[sp],g_size[sp]=parse_gff3.main(i,gene_size,gsize=True)
        except TypeError:
            pass
    # print(t_intro)
    return t_intro,g_size

def print_error(value):
    print("error: ", value)

def count_intro_and_gene(trans:pd,hgt_out:str,arg:dict):
    intro_base=arg['intro']
    gsize_base=arg['gsize']
    c_sp_list=arg['sp_list']
    tree_fp=arg['tree_fp']
    pattern = r"'([^'\s]*)'"
    HGT_out=open(hgt_out,'w')
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
        tree_file=os.path.join(tree_fp,og_id+"_ali.mafft_rmdup_trimal_JTT.treefile")
        gene_tree=Tree(tree_file)
        # tree_name=os.path.basename(tree_file)
        intro_out=open('./gene_stat_result/'+og_id+'_'+str(loc)+"_intro",'w')
        gsize_out=open('./gene_stat_result/'+og_id+'_'+str(loc)+"_gsize",'w')
        HGT_stat={'intro_num':[],'intro_size':[],'gene_size':[],'gene_loc':[]}
        no_sp=[]
        for i in gene_tree.get_leaf_names():
            sp=i.split('_')[0]
            if sp in c_sp_list:
                pass
            else:
                sp=[x for x in c_sp_list if x.startswith(sp[:4])]
                if not sp:
                    if i.split('_')[0] not in no_sp:
                        print(i.split('_')[0],'not include in gff3 file!')
                        no_sp.append(i.split('_')[0])
                    continue
                else:
                    sp=sp[0]
            if sp in intro_base:
                # print(i)
                c_intro_base=intro_base[sp]
                c_gsize_base=gsize_base[sp]
                c_gene_list=c_intro_base.keys()
                if i in c_gene_list:
                    tree_gene_name=i
                else:
                    try:
                        tree_gene_name=difflib.get_close_matches(i.split('_',maxsplit=1)[1],c_gene_list,1, cutoff=0.4)[0]
                    except IndexError:
                        print(i,'can\'t find proper genes')
                        continue
                #gene_id intro_num intro_size
                intro_out.write(i+'\t'+str(c_intro_base[tree_gene_name][0])+'\t'+str(c_intro_base[tree_gene_name][1]))
                intro_out.write('\n')
                #gene_id genesize gene_loc
                gsize_out.write(i+'\t'+str(c_gsize_base[tree_gene_name][0])+'\t'+str(c_gsize_base[tree_gene_name][1]))
                gsize_out.write('\n')
                if i in HGT_gene:
                    HGT_stat['intro_num'].append(c_intro_base[tree_gene_name][0])
                    HGT_stat['intro_size'].append(c_intro_base[tree_gene_name][1])
                    HGT_stat['gene_size'].append(c_gsize_base[tree_gene_name][0])
                    HGT_stat['gene_loc'].append(c_gsize_base[tree_gene_name][1])
        intro_out.close()
        gsize_out.close()
        HGT_out.write(str(loc)+"\t"+og_id+"\t"+",".join(HGT_gene)+"\t"+str(np.median(HGT_stat['intro_num']))\
            +"\t"+str(np.median(HGT_stat['intro_size']))+"\t"+str(np.median(HGT_stat['gene_size']))+"\t"+str(np.median(HGT_stat['gene_loc'])))
        HGT_out.write('\n')
        # print(datetime.now(),)

if __name__ == "__main__":
    #USAGE python trans_file.xlsx tree_fp gff3_file gene_size
    # trans_file=r"E:\Bio_analysis\Grass_horizatal_transversion\6_final_intergration\0724_manully_check.xlsx"
    trans_file=sys.argv[1]
    # tree_file_list=r"E:\Bio_analysis\Grass_horizatal_transversion\5_find_clade\OG0000805_ali.mafft_rmdup_trimal_JTT.treefile"
    tree_fp=sys.argv[2]
    # tree_fp=r'E:\Bio_analysis\Grass_horizatal_transversion\7_final_plot\0805manully_check_circular_tree_set'
    gff3_file=sys.argv[3]
    # gff3_file=r"E:\Bio_analysis\Grass_horizatal_transversion\8_make_additional\prase_gff\gff_name.list"
    # /public/home/huangyj/grass_horizatal/ref/ful_gff/Acom_full.gff3
    # /public/home/huangyj/grass_horizatal/ref/ful_gff/Aspl_full.gff3
    # /public/home/huangyj/grass_horizatal/ref/ful_gff/Astr_full.gff3
    # /public/home/huangyj/grass_horizatal/ref/ful_gff/Bamp_full.gff3
    # /public/home/huangyj/grass_horizatal/ref/ful_gff/Bdis_full.gff3
    gene_size=sys.argv[4]
    # gene_size=r"E:\Bio_analysis\Grass_horizatal_transversion\8_make_additional\genome.size"
    # sp1_chr1  100000
    # sp1_chr2  100000
    # sp1_chr3  100000
    # sp2_chr3  100000
    split_num=20
    # main(trans_file,tree_fp,gff3_file,gene_size,split_num=20)
    if not os.path.isdir('./gene_stat_result'):
        os.mkdir('./gene_stat_result')
    gff_file_name_list=[]
    print(11)
    with open(gff3_file,'r') as g_f:
        for i in g_f:
            line=i.strip()
            gff_file_name_list.append(line)
    # print(11)
    intro_base,gsize_base=get_gff_gene_base(gff_file_name_list,gene_size)
    # print(11)
    c_sp_list=[os.path.basename(x).split('_')[0] for x in gff_file_name_list]
    trans=pd.read_excel(trans_file,engine='openpyxl',sheet_name=1)

    p=multiprocessing.Pool()
    multiprocessing.log_to_stderr()
    arg_dic={'intro':intro_base,'gsize':gsize_base,'sp_list':c_sp_list,'tree_fp':tree_fp}
    print(datetime.now(),'Parent process %s.' % os.getpid())
    for x in range(0,split_num):
        c_index=[y for y in trans.index.tolist() if int(y) % int(split_num) == x]
        c_dataframe=trans.iloc[c_index]
        c_hgt="HGT_subprocesses"+str(x)+'.txt'
        p.apply_async(LogExceptions(count_intro_and_gene), args=(c_dataframe,c_hgt,arg_dic, ))
    print(datetime.now(),'Waiting for all subprocesses done...')
    p.close()
    p.join()
    print(datetime.now(),'All subprocesses done!')