'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-10-10 20:03:13
LastEditors: Ne0tea
LastEditTime: 2024-12-13 16:45:08
'''
from ete3 import Tree
import sys
import src.para_base  as pb
import src.global_variable as glv
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import multiprocessing
from scipy import stats
def define_sub(gene,clade_subg_dic):
    sp=gene.split("_")[0]
    for i in clade_subg_dic:
        if sp in clade_subg_dic[i]:
            return i

def get_subfamily_list(tree:Tree):
    root = tree
    if len(root.children) == 2:
        # 获取两个分组
        group1 = root.children[0]
        group2 = root.children[1]

        # 提取两个分组的叶子节点
        subfamily1 = [leaf.name for leaf in group1.iter_leaves()]
        subfamily2 = [leaf.name for leaf in group2.iter_leaves()]
    return [subfamily1, subfamily2]

def get_subfamily_list_from_phy(temp_phylo,OUT):
    if not OUT: #if there is no OUT in sub phylogeny, just return tree
        pass
    else:
        try:
            new_root = temp_phylo&OUT
            temp_phylo.set_outgroup(new_root)
            sister_root = new_root.get_sisters()[0]
            node=temp_phylo&OUT
            node.delete()
            temp_phylo=sister_root
        except:
            return get_subfamily_list(temp_phylo)
    return get_subfamily_list(temp_phylo)

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

def get_receive_donor_distances(gf_id, og_id, gf_type,tree_nwk:str, receive_gene_list, donor_gene_list, donor, receive, subfamily_list, clade_subg_dic):

    tree = Tree(tree_nwk)
    target_gene_list = receive_gene_list.split('|')
    receive_sp_list = list(set([i.split('_')[0] for i in target_gene_list]))
    donor_gene_list = donor_gene_list.split('|')
    same_subf = subfamily_list[0] if receive in subfamily_list[0] else subfamily_list[1]
    clade_gene_dic = {}
    native_receive_sp = []
    native_copy = 'False'
    for i in tree.get_leaf_names():
        gene_clade=define_sub(i,clade_subg_dic)
        if i in target_gene_list or i in donor_gene_list : continue
        if gene_clade not in subfamily_list[0] and gene_clade not in subfamily_list[1]: continue
        if gene_clade == receive:
            native_receive_sp.append(i.split('_')[0])
        if gene_clade in clade_gene_dic:
            clade_gene_dic[gene_clade].append(i)
        else:
            clade_gene_dic[gene_clade]=[i]

    if len(set(receive_sp_list) & set(native_receive_sp)) >= 3:
        native_copy = 'True'

    distances = {}
    for rg in target_gene_list:
        rg_node = tree&rg
        every_donor_list=[]
        for dg in donor_gene_list:
            if dg in tree.get_leaf_names():
                dg_node = tree&dg
            else:
                continue
            cur_distance = tree.get_distance(rg_node,dg_node)
            every_donor_list.append(cur_distance)
            if 'Receive' in distances:
                distances['Receive'].append(cur_distance)
            else:
                distances['Receive'] = [cur_distance]
    for dg in donor_gene_list:
        if dg in tree.get_leaf_names():
            dg_node = tree&dg
        else:
            continue
        every_donor_list=[]
        for clade in same_subf:
            if clade not in clade_gene_dic: continue
            # print(clade)
            for sg in clade_gene_dic[clade]:
                sg_node = tree&sg
                cur_distance = tree.get_distance(sg_node,dg_node)
                every_donor_list.append(cur_distance)
                if 'Sister' in distances:
                    distances['Sister'].append(cur_distance)
                else:
                    distances['Sister'] = [cur_distance]

    return gf_id, og_id, gf_type, distances, native_copy

def plot_distance_distribution(distances,template_color,pre='distance_freq'):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'

    plt.figure(figsize=(6, 4))

    for cur_clade in distances:
        plt.hist(distances[cur_clade], bins=75, histtype='step', label=cur_clade,edgecolor=template_color[cur_clade])
    plt.legend()
    plt.title('')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.savefig(pre+'.pdf')
    # plt.show()
    plt.close()

def run_test_batch(opt_list:list,config, num_processes, outfile):
    print(f'Modified branch length test running start with {len(opt_list)}')
    pb.main(config)
    sub_phlo=glv.get('sub_phlo')
    clade_subg_dic=glv.get('clade_subg_dic')
    out_name = decide_out_name(clade_subg_dic)
    subfamily_list=get_subfamily_list_from_phy(Tree(sub_phlo),out_name)

    p_values = []
    out_line=[]
    out_put_list = []
    ##(gf_id, og_id, gf_type, cur_node,gene_list,donor_gene_list,donor,receiver)
    para_pass_multi=[(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],subfamily_list, clade_subg_dic) for x in opt_list]

    with multiprocessing.Pool(processes=num_processes) as pool:
        gf_id_list, og_id_list, gf_type_list, distance_dic_list, native_copy_list = [], [], [], [], []
        parallel_result = pool.starmap(get_receive_donor_distances, para_pass_multi)
        for i in parallel_result:
            gf_id_list.append(i[0])
            og_id_list.append(i[1])
            gf_type_list.append(i[2])
            distance_dic_list.append(i[3])
            native_copy_list.append(i[4])

    for idx, distance_dic in enumerate(distance_dic_list):
        native_copy = native_copy_list[idx]
        if 'Sister' not in distance_dic:
            output = '\t'.join([gf_id_list[idx],og_id_list[idx],gf_type_list[idx], 'False','NA','NA','NA']) + '\tNA'
            out_put_list.append(output)
            continue
        stat, p_value = stats.ttest_ind(distance_dic['Receive'], distance_dic['Sister'], permutations=100,alternative='less', nan_policy='omit')
		# stat, p_value = stats.mannwhitneyu(distance_dic['Receive'], distance_dic['Sister'], permutations=100,alternative='less', nan_policy='omit')
        ## effect size calculation
        n1 = len(distance_dic['Receive'])
        n2 = len(distance_dic['Sister'])
        z_value = stats.norm.ppf(p_value / 2)
        effect_size = z_value / ((n1 + n2) ** (1/2))
        p_values.append(p_value)
        out_line.append([gf_id_list[idx],og_id_list[idx],gf_type_list[idx],native_copy,str(z_value),str(effect_size),str(p_value)])

    if len(p_values) == 0:
        print('No available p-value, Skip FDR. Return.')
        return

    _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
    for out_unit, q_value in zip(out_line, q_values):
        output = '\t'.join(out_unit) + '\t' + str(q_value)
        out_put_list.append(output)
    outf = open(outfile, 'w')
    for i in out_put_list:
        outf.write(i+'\n')
    outf.close()

if __name__ == "__main__":
    agfd_outfile=sys.argv[1]
    # outfile = r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\test.out'
    #GF_8    OG0000005       1       RGF     PAN     ORY     Pvir-N|Cmac-A|Cmac-B    Oalt-D  Oalt-D_prot_Oalt-D_OalD10g102020_1
    # nwk_file = r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\test.out.OG.node.nwk'
    nwk_file = sys.argv[2]
    # config = r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt'
    config=sys.argv[3]
    outputfile=sys.argv[4]
    template_color={"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1",'Receive':'lightgrey','Sister':'black'}

    nwk_dic = {}
    with open(nwk_file, 'r') as nf:
        for line in nf:
            line = line.strip().split()
            nwk_dic[line[0]] = line[3]

    opt_list=[]
    with open(agfd_outfile, 'r') as of:
        for line in of:
            line=line.strip()
            gf_id, og_id, gene_len, gf_type, donor, receiver, _, donor_gene_list, _, gene_list = line.split('\t')
            opt_list.append((gf_id, og_id, gf_type, nwk_dic[gf_id], gene_list, donor_gene_list, donor, receiver))
    run_test_batch(opt_list, config, 10, outputfile)