'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-29 16:56:20
LastEditors: Ne0tea
LastEditTime: 2024-11-13 00:38:11
'''
import argparse
import pandas as pd
import time

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

def main(outfile, config_file):
    clade_subg_dic={}
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
    sp_ord=[]
    for x in clade_subg_dic:
        if x != decide_out_name(clade_subg_dic):
            sp_ord.extend(clade_subg_dic[x])
    clade_ord=[ x for x in clade_subg_dic if x != decide_out_name(clade_subg_dic)]
    Donor_receive_spscale_df=pd.DataFrame(index=sp_ord,columns=sp_ord)
    Donor_receive_Cladescale_df=pd.DataFrame(index=clade_ord,columns=clade_ord)

    with open(outfile,'r') as outf:
        for line in outf:
            cur_line=line.strip().split('\t')
            cur_donor_clade, cur_receive_clade, cur_donor_sp, donor_gene, cur_receive_sp=cur_line[4:9]
            single_donor_sp=cur_donor_sp.split('|')
            single_receive_sp=cur_receive_sp.split('|')
            for receive in single_receive_sp:
                if receive not in sp_ord:
                    continue
                for donor in single_donor_sp:
                    if donor not in sp_ord:
                        continue
                    if not pd.isna(Donor_receive_spscale_df.loc[receive,donor]):
                        Donor_receive_spscale_df.loc[receive,donor]+=1
                    else:
                        Donor_receive_spscale_df.loc[receive,donor]=1

            single_donor_clade=cur_donor_clade.split('|')
            single_receive_clade=cur_receive_clade.split('|')

            for receive in single_receive_clade:
                if receive not in clade_ord:
                    continue
                for donor in single_donor_clade:
                    if donor not in clade_ord:
                        continue
                    if not pd.isna(Donor_receive_Cladescale_df.loc[receive,donor]):
                        Donor_receive_Cladescale_df.loc[receive,donor]+=1
                    else:
                        Donor_receive_Cladescale_df.loc[receive,donor]=1
    Donor_receive_Cladescale_df.to_csv(outfile+'_clade.csv',index=True,header=True)
    Donor_receive_spscale_df.to_csv(outfile+'_sp.csv',index=True,header=True)

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
    GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
    '''
    outfile=r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out.consistance.significant.all.out'
    config_file=r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    main(outfile, config_file)

    end = time.perf_counter()
    h=int((end - begin)/3600)
    m=int((end - begin)/60 - h*60)
    s=int(end - begin - h*3600 - m*60)
    print("All tasks used time: %sh %sm %ss" % (h,m,s))