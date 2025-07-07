'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-10-03 15:17:05
LastEditors: Ne0tea
LastEditTime: 2024-11-13 00:41:46
'''
'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-29 16:56:20
LastEditors: Ne0tea
LastEditTime: 2024-10-16 10:25:18
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

def main(outfile):
    receive_from_phy_stat_output=outfile+'.receive_from_phy_stat'
    donor_to_phy_stat_output=outfile+'.donor_to_phy_stat'
    receive_from_donor_dic={}
    donor_to_receive_dic={}
    '''
    out file format example:
    donor_clade, receive_clade, donor_sp, receive_sp
    GF_6    OG0000005       30      Toplogy inconsistent    BAM     CHL     Dlat-B|Pedu-C...       Aluo-D_prot_Aluo-D_Alu06Dg13070_mRNA1|Pedu-C_prot_Pedu-C_EVM0028367_1...    Etef-B|Cson-B...  Ecur_prot_Ecur_TVU20393_m|Cson-B_prot_Cson-B_CsB501278_1...
    '''
    with open(outfile,'r') as outf:
        for i in outf:
            cur_line=i.strip().split('\t')
            donor_clade, receive_clade, donor_sp, donor_gene, receive_sp=cur_line[4:9]
            receive_sp_list=list(set(receive_sp.split('|')))
            receive_sp_list.sort()
            donor_sp_list=list(set(donor_sp.split('|')))
            donor_sp_list.sort()
            receive_sp_join='|'.join(receive_sp_list)
            donor_sp_join = '|'.join(donor_sp_list)
            donor_clade_sep = donor_clade.split('|')
            receive_clade_sep = receive_clade.split('|')
            for single_donor in donor_clade_sep:
                if receive_sp_join in receive_from_donor_dic:
                    if single_donor in receive_from_donor_dic[receive_sp_join]:
                        receive_from_donor_dic[receive_sp_join][single_donor]+=1
                    else:
                        receive_from_donor_dic[receive_sp_join][single_donor]=1
                else:
                    receive_from_donor_dic[receive_sp_join]={single_donor:1}
            for single_receive in receive_clade_sep:
                if donor_sp_join in donor_to_receive_dic:
                    if single_receive in donor_to_receive_dic[donor_sp_join]:
                        donor_to_receive_dic[donor_sp_join][single_receive]+=1
                    else:
                        donor_to_receive_dic[donor_sp_join][single_receive]=1
                else:
                    donor_to_receive_dic[donor_sp_join]={single_receive:1}

    sp_receive_from_donor_df = pd.DataFrame.from_dict(receive_from_donor_dic, orient='index').reset_index()
    sp_donor_to_receive_df = pd.DataFrame.from_dict(donor_to_receive_dic, orient='index').reset_index()
    sp_receive_from_donor_df.fillna(0, inplace=True)
    sp_receive_from_donor_df.rename(columns={'index': 'Phylogeny'}, inplace=True)
    sp_donor_to_receive_df.fillna(0, inplace=True)
    sp_donor_to_receive_df.rename(columns={'index': 'Phylogeny'}, inplace=True)

    sp_receive_from_donor_df.to_csv(receive_from_phy_stat_output,index=False,sep='\t')
    sp_donor_to_receive_df.to_csv(donor_to_phy_stat_output,index=False,sep='\t')


if __name__ == "__main__":
    begin = time.perf_counter()
    # parser = argparse.ArgumentParser(description="Stat donor/receiver in phylogeny format from the AGFD result.")
    # parser.add_argument("-i","--outfile", type=str, help="Outfile from AGFD.")
    # parser.add_argument("-c","--config", type=str, help="Path to the config file.")
    # args = parser.parse_args()
    # main(args.outfile, args.config)
    '''
    out file format example:
    donor_clade, receive_clade, donor_sp, receive_sp
    OG0000020	Sbic_prot_Sbic_Sobic_004G235600_1_v3_2|Eoph_prot_Eoph_ctg632_38_m	2   ORY	PAN	Ogla|Oruf|Osat	Eoph|Sbic
    '''
    outfile=r'E:\Bio_analysis\HGT_newcollection\4_transfer_destiny\Poaceae_tree_V16.out.consistance.significant.all.out'
    main(outfile)
    end = time.perf_counter()
    h=int((end - begin)/3600)
    m=int((end - begin)/60 - h*60)
    s=int(end - begin - h*3600 - m*60)
    print("All tasks used time: %sh %sm %ss" % (h,m,s))