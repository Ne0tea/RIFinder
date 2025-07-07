'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-10-20 17:18:40
LastEditors: Ne0tea
LastEditTime: 2024-10-15 14:06:05
'''
import src.global_variable as glv

def main(config:str):
    glv._init_()
    clade_subg_dic={}

    with open(config,'r') as c_file:
        phlo=0
        subg=0
        sp=0
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


    glv.set('clade_subg_dic',clade_subg_dic)
    glv.set('sub_phlo',sub_phlo)
    # glv.set('gene_order',gene_order)

if __name__ == "__main__":
    config=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    gene_order=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\All_gene.ord'
    main(config,gene_order)