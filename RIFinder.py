#!/bin/env python
'''
Descripttion: RIFinder main progress
Author: Ne0tea
version: V0.1.0
Date: 2023-07-18 23:54:12
LastEditors: Ne0tea
LastEditTime: 2025-02-15 19:44:23
'''
import re
import argparse
import time
import src.major_detect as md
import traceback
from src.log import get_logger, close_log_queue
import warnings
import os
import multiprocessing
import src.modified_branch_lenth_test as mblt
import src.para_base as pb
import src.global_variable as glv
import src.get_prop_donor as get_prop_donor
def run_command(gene_tree,sup_value,config,rootnested):

    pb.main(config)
    Dect_result=[]
    top_node_result=[]
    if rootnested:
        Ngenenode_set_result=[]
    try:
        id=re.search('OG[0-9]{7}',gene_tree).group()
    except AttributeError:
        id=os.path.basename(gene_tree)
        logger.info(id+" isn't built by orthofinder. ")
    logger.info('tree id {0}, support value {1}'.format(id,sup_value))
    logger.info(id + ' analysis start')
    if rootnested:
        a_HGT_gene,a_HGT_gene_len,a_HGT_topnode_dic,clade_state,Ngenenode_set, r_GF_dic=md.main(gene_tree,sup_value,rootnested)
        Ngenenode_set_result.append(Ngenenode_set)
    else:
        a_HGT_gene,a_HGT_gene_len,a_HGT_topnode_dic,clade_state, r_GF_dic=md.main(gene_tree,sup_value,rootnested)

    if clade_state:
        # logger.info("{0} identified ancient gene flow genes".format(id))
        logger.info(id + " identified ancient gene flow genes")
        for i in a_HGT_gene:
            status, donor_clade, receive_clade, donor_set, receive_sp=get_prop_donor.main(a_HGT_gene[i],a_HGT_topnode_dic[i])
            if donor_set[0]:
                out_line="\t".join([id, str(a_HGT_gene_len[i]), status, donor_clade, receive_clade, donor_set[0], donor_set[1], receive_sp, "|".join(a_HGT_gene[i])])
                Dect_result.append(out_line)
                node_out_line="\t".join([id,str(a_HGT_gene_len[i]),a_HGT_topnode_dic[i].write()])
                top_node_result.append(node_out_line)
    else:
        logger.info(id+" detected no ancient gene flow")

    for rrecord in r_GF_dic:
        donor_clade, donor_sp_set, receive_clade, receive_sp, donor_gene, receive_gene, sg_newick = r_GF_dic[rrecord]
        out_line="\t".join([id, str(len(receive_gene)), "RGF", donor_clade, receive_clade, donor_sp_set, donor_gene, receive_sp, "|".join(receive_gene)])
        Dect_result.append(out_line)
        node_out_line="\t".join([id, str(len(receive_gene)), sg_newick])
        top_node_result.append(node_out_line)

    logger.info(id+' analysis Done')
    if rootnested:
        return Dect_result, top_node_result, Ngenenode_set_result
    else:
        return Dect_result, top_node_result

def main(args):
    config=args.config

    gene_tree_set_file = args.tree
    sup_value=args.support
    num_processes=args.thread
    rootnested=args.rootnested
    out=args.out

    global logger
    logger = get_logger(out+'.log','RIFinder')

    error_flag = multiprocessing.Value('i', 0)

    pb.main(config)
    clade_subg_dic=glv.get("clade_subg_dic")
    print('Porcess with clade as {0}'.format(list(clade_subg_dic.keys())))

    with open(gene_tree_set_file,'r') as g_file:
        parallel_command_value=[(x.strip(),sup_value,config,rootnested) for x in g_file]

    try:
        with multiprocessing.Pool(processes=num_processes) as pool:
            if rootnested:
                AGF_result, AGF_top_node_result, AGF_Ngenenode_result = [], [], []
                parallel_result = pool.starmap(run_command, parallel_command_value)
                for i in parallel_result:
                    AGF_result.extend(i[0])
                    AGF_top_node_result.extend(i[1])
                    AGF_Ngenenode_result.extend(i[2])

            else:
                AGF_result, AGF_top_node_result = [], []
                parallel_result = pool.starmap(run_command, parallel_command_value)
                for i in parallel_result:
                    AGF_result.extend(i[0])
                    AGF_top_node_result.extend(i[1])

    except Exception as e:
        tb_info = traceback.format_exc()
        print(f"An error occurred: {e}")
        print(f"Traceback information: {tb_info}")
        error_flag.value = 1

    if error_flag.value:
        # log_queue.put(None)
        # listener_process.join()
        raise RuntimeError("Failed to execute tasks.")
    else:
        print(f"All task completed. {len(parallel_command_value)} gene trees were disentangled.")
        out_file=open(out,'w')
        top_node_out_file=open(out+'.OG.node.nwk','w')
        if rootnested:
            ngout_file=open('Ngenenode_set.txt','w')

        GF_count=0
        for jobs in AGF_result:
            GF_count+=1
            out_file.write('GF_'+str(GF_count)+'\t')
            out_file.write(jobs+'\n')

        GF_count=0
        for jobs in AGF_top_node_result:
            GF_count+=1
            top_node_out_file.write('GF_'+str(GF_count)+'\t')
            top_node_out_file.write(jobs+'\n')

        if rootnested:
            for jobs in AGF_Ngenenode_result:
                for result in jobs:
                    for line in result:
                        ngout_file.write(line+'\n')
            ngout_file.close()

        out_file.close()
        top_node_out_file.close()

        ###Make introgression pvalue test
        p_value_outfile = out+'.p_value'
        nwk_dic = {}
        with open(out+'.OG.node.nwk', 'r') as nf:
            for line in nf:
                line = line.strip().split()
                nwk_dic[line[0]] = line[3]

        opt_list=[]
        with open(out, 'r') as of:
            for line in of:
                line=line.strip()
                if line.split('\t')[3] != 'RGF':
                    gf_id, og_id, gene_len, gf_type, donor, receiver, _, donor_gene_list, _, gene_list = line.split('\t')
                    opt_list.append((gf_id,og_id,gf_type,nwk_dic[gf_id],gene_list,donor_gene_list,donor,receiver))
        if opt_list:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")  # get all warnings
                mblt.run_test_batch(opt_list, config, num_processes, p_value_outfile)
                if w:
                    for warning in w:
                        logger.warning(f"Warning: {warning.message}")
    close_log_queue()

if __name__ == "__main__":
    begin = time.perf_counter()
    parser = argparse.ArgumentParser(prog="RIFinder",formatter_class=argparse.RawTextHelpFormatter,description="""
-------------------------------------------------------------------------------------------------------
RIFinder
Author: Yujie Huang <yujiehuang@zju.edu.cn>, ZJU
        Dongya Wu wudongya@zju.edu.cn,ZJU
Version: v0.0.1
Identify remote introgression across a phylogeny of multiple assemblies.
-------------------------------------------------------------------------------------------------------
""")
    parser.add_argument('-f','--treefile', dest='tree', help="a NWK gene tree file set\n##one tree file per line", type=str,required=True)
    parser.add_argument('-c','--config',help="a config file indicating phlogentic relationship and subgroup species\n##Please check example config file for exact content", 
                        type=str,required=True)
    parser.add_argument('-o','--out',help="output alignments to FILE [RIFinder_output.txt]", type=str,default='RIFinder_output.txt')
    parser.add_argument('-rn','--rootnested',help="Whether or not to disentangle tree structure nested with Outgroup [Flase]", action='store_true')
    parser.add_argument('-t',dest='thread',type=int,default=1,help="number of threads [1] ")
    parser.add_argument('--support',type=int,default=70,help="Filter out branches with bootstrap values lower than the standard during the identification of ancient transfer events. [70/0.7] ")
    args = parser.parse_args()

    main(args)
    end = time.perf_counter()
    h=int((end - begin)/3600)
    m=int((end - begin)/60 - h*60)
    s=int(end - begin - h*3600 - m*60)
    print("All tasks used time: %sh %sm %ss" % (h,m,s))