'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-10-20 17:18:40
LastEditors: Ne0tea
LastEditTime: 2023-10-20 22:42:16
'''
import global_variable as glv

def main(config:str,gene_order:str):
    glv._init()
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
    glv.set('gene_order',gene_order)

if __name__ == "__main__":
    main()