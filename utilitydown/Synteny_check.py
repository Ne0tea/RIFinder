'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-10-15 16:17:13
LastEditors: Ne0tea
LastEditTime: 2024-10-15 20:36:15
'''
import time
import re

def read_config(config_file):
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
    return clade_subg_dic,subg_clade_dic

def main(raw_out_file, config_file, syntenic_gene_pan_file):
    clade_subg_dic, subg_clade_dic= read_config(config_file)
    outfile = open(raw_out_file+'nosyn.out','w')
    all_gene_sg = {}
    with open(syntenic_gene_pan_file,'r') as synf:
        for i in synf:
            cur_line=i.strip().split()[4:]
            cur_sg = i.strip().split()[0]
            for genes in cur_line:
                for gene in genes.split(','):
                    cur_gene_clade = subg_clade_dic[gene.split('_')[0]]
                    sp = gene.split('_')[0]
                    gene = gene.split('_', 1)[-1]
                    all_gene_sg[gene] = cur_sg
    with open(raw_out_file,'r') as outf:
        for line in outf:
            cur_line=line.strip().split('\t')
            '''
            out file format example:
            donor_clade, receive_clade, donor_sp, receive_sp
            GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
            '''
            gf_id, donor_clade,donor_gene,receive_gene=cur_line[0], cur_line[4],cur_line[7],cur_line[9]
            receive_sg_list = []
            for gene in receive_gene.split('|'):
                sp_name = gene.split('_')[0]
                gene_name = re.search(r'_prot_(.*)', gene).group(1)
                if gene_name in all_gene_sg:
                    cur_sg = all_gene_sg[gene_name]
                else:
                    cur_sg = 'Left'
                receive_sg_list.append(cur_sg)
            donor_sg_list = []
            for gene in donor_gene.split('|'):
                sp_name = gene.split('_')[0]
                gene_name = re.search(r'_prot_(.*)', gene).group(1)
                if gene_name in all_gene_sg:
                    cur_sg = all_gene_sg[gene_name]
                else:
                    cur_sg = 'Left'
                donor_sg_list.append(cur_sg)
            intersection = set(receive_sg_list) & set(donor_sg_list)
            if len(intersection) > 1 and ('left' not in intersection):
                print(gf_id)
            else:
                outfile.write(line)

if __name__ == "__main__":
    begin = time.perf_counter()
    # parser = argparse.ArgumentParser(description="Stat donor/receiver from the AGFD result.")
    # parser.add_argument("-i","--outfile", type=str, help="Outfile from AGFD.")
    # parser.add_argument("-c","--config", type=str, help="Path to the config file.")
    # args = parser.parse_args()
    # main(args.outfile, args.config)
    '''
    out file format example:
    donor_clade, receive_clade, donor_sp, receive_sp
    OG0000020	Sbic_prot_Sbic_Sobic_004G235600_1_v3_2|Eoph_prot_Eoph_ctg632_38_m	2   ORY	PAN	Ogla|Oruf|Osat	Eoph|Sbic
    '''

    raw_out_file=r'E:\Bio_analysis\HGT_newcollection\3_outfile_stat\Poaceae_tree_V14.out.consistance.significant.out'
    config_file=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    syntenic_gene_pan_file = r'E:\Bio_analysis\HGT_newcollection\3_outfile_stat\120_sub_m108.list.SG.pan'

    main(raw_out_file, config_file, syntenic_gene_pan_file)

    end = time.perf_counter()
    h=int((end - begin)/3600)
    m=int((end - begin)/60 - h*60)
    s=int(end - begin - h*3600 - m*60)
    print("All tasks used time: %sh %sm %ss" % (h,m,s))