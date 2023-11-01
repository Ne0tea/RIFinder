'''
Descripttion: HGTfinder main progress
Author: Ne0tea
version: V1.0.0
Date: 2023-07-18 23:54:12
LastEditors: Ne0tea
LastEditTime: 2023-10-21 13:59:50
'''
import re
import argparse
import time
import collapse_filter
import structure_extract_further
import find_clade_v7
import tree_KNN_cluster_v3
import logging
import os
import para_base  as pb
import global_variable as glv

# def run_command(cmd):
#     # print("INFO: Running command: {0}".format(cmd), flush=True)
#     print(cmd)
#     return_code = subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL)
#     if return_code != 0:
#         print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))


def main(args):
    gene_tree_set_file = args.tree
    dom_value=args.kmax
    minor_value=args.kmin
    sup_value=args.support
    out=args.out
    out_file=open(out,'w')
    with open(gene_tree_set_file,'r') as g_file:
        gene_tree_set=[x.strip() for x in g_file]
    for gene_tree in gene_tree_set:
        try:
            og=re.search('OG[0-9]{7}',gene_tree).group()
            logging.info('tree id {0}, dominate value {1}, minor value {2}, support value {3}'.format(og,dom_value,minor_value,sup_value))
            id=og
        except AttributeError:
            print("WARING: gene tree name can't be recongnized !")
            logging.info('tree id {0}, dominate value {1}, minor value {2}, support value {3}'.format('tree1',dom_value,minor_value,sup_value))
            id='tree1'
        logging.info('Start {0} analysis'.format(id))
        print(id,'analysis start')
        logging.info("Step1: filter gene tree with inconsistent structures.")
        inc_reason,inc_state=collapse_filter.main(gene_tree)
        logging.info("-------------------------------------------------------------------------------------------------------")
        logging.info("Step2: filter gene tree with mutiple collapsed leaves.")

        collapsed_leaf_num,col_state=structure_extract_further.main(gene_tree)
        if inc_state and col_state:
            r_HGT=True
            logging.info("-------------------------------------------------------------------------------------------------------")
            logging.info('Success: {0} may be recent HGT'.format(id))
            print(id,'may be recent HGT')
        elif not inc_state:
            r_HGT=False
            logging.info("-------------------------------------------------------------------------------------------------------")
            logging.info('WARING: {0} has too much inconsistent gene as {1}'.format(id,inc_reason))
            logging.info("Skip Step3: recent HGT identification")
            print("Skip Step3: recent HGT identification")
        elif not col_state:
            r_HGT=False
            logging.info("-------------------------------------------------------------------------------------------------------")
            logging.info('WARING: {0} has too much collapsed leaves {1}'.format(id,collapsed_leaf_num))
            logging.info("Skip Step3: recent HGT identification")
            print("Skip Step3: recent HGT identification")

        if r_HGT:
            logging.info("-------------------------------------------------------------------------------------------------------")
            logging.info("Step3: identify recent HGT gene in the tree.")
            r_HGT_gene,KNN_state=tree_KNN_cluster_v3.main(gene_tree,dom_value,minor_value,sup_value)
            if KNN_state:
                logging.info("{0} identified recent HGT genes".format(id))
                print(id,"identified recent HGT genes")
                for i in r_HGT_gene:
                    out_line="\t".join(["G_line_R",id,str(r_HGT_gene[i]),str(len(r_HGT_gene[i]))])
                    out_file.write(out_line+'\n')
                    # print("G_line_R",og,r_HGT_gene[i],len(r_HGT_gene[i]),sep='\t')
            else:
                print(id,"identified no recent HGT gene")
        logging.info("-------------------------------------------------------------------------------------------------------")
        logging.info("Step4: identify ancient HGT gene in the tree.")
        a_HGT_gene,a_HGT_gene_len,clade_state=find_clade_v7.main(gene_tree,sup_value)
        # a_HGT_gene,a_HGT_gene_len,clade_state=main(gene_tree,sup_value)
        if clade_state:
            logging.info("{0} identified ancient HGT genes".format(id))
            print(id,"identified ancient HGT genes")
            for i in a_HGT_gene:
                out_line="\t".join(["G_line_A",id,str(a_HGT_gene[i]),str(a_HGT_gene_len[i])])
                out_file.write(out_line+'\n')
                # print("G_line_A",og,a_HGT_gene[i],a_HGT_gene_len[i],sep='\t')
        else:
            print(id,"identified no ancient HGT gene")
        logging.info('{0} analysis done'.format(id))
        logging.info('')
        print(id,'analysis Done')

if __name__ == "__main__":
    begin = time.perf_counter()
    parser = argparse.ArgumentParser(prog="HgtFinder",formatter_class=argparse.RawTextHelpFormatter,description="""
-------------------------------------------------------------------------------------------------------
HgtFinder
Author: Yujie Huang <12116008@zju.edu.cn>, ZJU
Version: v1.0
Identify HGT events based on gene tree
-------------------------------------------------------------------------------------------------------
""")
    parser.add_argument('-t','--tree', help="supply a NWK gene tree file set", type=str,required=True)
    parser.add_argument('-c','--config',help="supply config file indicating phlogenic relationship and clade name", type=str,required=True)
    parser.add_argument('-go','--gene-order',help="supply gene order file", dest='go',type=str,required=True)
    parser.add_argument('-o','--out',help="supply out info file", type=str,required=True)
    parser.add_argument('--kmax',type=float,default=0.8,help="the porpotion of dominate family in KNN cluster")
    parser.add_argument('--kmin',type=float,default=0.1,help="the porpotion of minor family in KNN cluster")
    parser.add_argument('--support',type=int,default=50,help="Filter out branches with bootstrap values lower \
than the standard during the identification of ancient transfer events.")
    args = parser.parse_args()
    gene_order=args.go
    config=args.config
    pb.main(config,gene_order)
    clade_subg_dic=glv.get("clade_subg_dic")

    logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s %(levelname)s %(message)s',
                    datefmt='%a %d %b %Y %H:%M:%S',
                    filename='log_HGTfinder',
                    filemode='w')

    print('Porcess with clade as {0}'.format(list(clade_subg_dic.keys())))
    main(args)
    end = time.perf_counter()
    print("All tasks used time: %ss" % (end - begin))
