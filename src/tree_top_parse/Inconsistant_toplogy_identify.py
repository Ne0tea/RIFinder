'''
Descripttion: Determining the events of OG-like structures
Author: Ne0tea
version: 
Date: 2024-09-12 10:01:15
LastEditors: Ne0tea
LastEditTime: 2024-12-13 14:13:40
'''

def decide_out_name(clade_subg_dic):
    if 'OUT' in clade_subg_dic:
        return 'OUT'
    elif 'Out' in clade_subg_dic:
        return 'Out'
    elif 'Outgroup' in clade_subg_dic:
        return 'Outgroup'
    elif 'outgroup' in clade_subg_dic:
        return 'outgroup'

'''
contary_step: indicating hierarchical level
TGF: Trickle gene flow
AGF: Ancient gene flow
'''
def rooted_toplogy(process,subfam_list,out_name:str,cur_clade_process):
    subf1,subf2=subfam_list
    contary_switch=False
    contary_count=0
    contary_step=0
    last_value=[]
    plus_not_None=0
    for index,value in enumerate(process):
        plus_value=value[len(last_value):]
        if not plus_value: # Skip empty plus
            continue
        plus_not_None+=1
        if plus_not_None==2:
            return_index=index-1
            for back_clade in cur_clade_process[:index-1][::-1]:
                if back_clade != cur_clade_process[0]:
                    return_index-=1
                else:
                    break

        cur_process_clade=['subf1' if x in subf1 else 'subf2' for x in plus_value if x != out_name]
        last_value=value
        if plus_not_None == 1:
            cur_clade='subf1' if plus_value[0] in subf1 else 'subf2'
            contary_clade='subf2' if plus_value[0] in subf1 else 'subf1'
            contary_clade_count = len(subf1) if contary_clade == 'subf1' else len(subf2)

            #Use Cayley for calculate candiate top number
            contary_clade_top_count=contary_clade_count**(contary_clade_count-2)
            continue
        if contary_clade in cur_process_clade and not contary_switch:
            if plus_not_None != 2:
                contary_switch = True
                contary_step += contary_clade_top_count / cur_process_clade.count(contary_clade) / len(cur_process_clade)
                if cur_clade in cur_process_clade:
                    break
            else:
                if cur_clade not in cur_process_clade:
                    contary_step += 1 / 2 / cur_process_clade.count(contary_clade) / len(cur_process_clade)
            contary_count+=cur_process_clade.count(contary_clade)
            continue

        if contary_switch:
            if contary_clade in cur_process_clade:
                contary_step += contary_clade_top_count / cur_process_clade.count(contary_clade) / len(cur_process_clade)
            if cur_clade in cur_process_clade:
                break
            else:
                contary_count+=cur_process_clade.count(contary_clade)

    # print(contary_count, contary_step, contary_clade_top_count)
    if contary_count >= contary_clade_count:
        if contary_step < contary_clade_top_count:
            return 'TGF', cur_clade, return_index
        else:
            return 'AGF', contary_clade, return_index
    else:
        return 'Inconstant_clade', None, None

def decide_toplogy_from_process(cur_ample_process,subfamily_list,out_name,cur_clade_process):
    first_process= [xx for xx in cur_ample_process if xx]
    if not first_process:
        return 'Outlier_clade', None, None
    first_process = first_process[0]
    if len(first_process) != 1:
        return 'Outlier_clade', None, None
    # print(cur_ample_process)
    if out_name in cur_ample_process[-1]:
        toplogy,clade_name,clade_index = rooted_toplogy(cur_ample_process,subfamily_list,out_name,cur_clade_process)
        # print(clade_index)
    else:
        # toplogy,clade_name,clade_index=unroot_toplogy(cur_ample_process,subfamily_list,out_name)
        return 'Norooted', None, None
    return toplogy,clade_name,clade_index
