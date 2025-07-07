'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-05-14 21:56:06
LastEditors: Ne0tea
LastEditTime: 2024-10-28 21:08:35
'''
from ete3 import Tree
import random
import argparse

template_color={"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1","OUT":"#a9a29c"}
def generate_random_colors(n):
    colors = []
    for _ in range(n):
        # 随机生成 RGB 颜色值
        color = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        colors.append('#{:02X}{:02X}{:02X}'.format(*color))
    return colors

def define_sub(gene,subg_clade_dic):
    sp=gene.split("_")[0]
    return subg_clade_dic[sp]

def write_itol_anno(treefile,color_template,subg_clade_dic):
    anno_out=open(treefile+"_itol_annotation.txt",'w')

    anno_out.write("DATASET_COLORSTRIP\nSEPARATOR\tTAB\n")
    anno_out.write("BORDER_WIDTH\t0.5\n")
    anno_out.write("COLOR\t#bebada\n")
    anno_out.write("DATASET_LABEL\tSub_family_divident\n")

    color_line="\t".join(color_template.values())
    label_line="\t".join(color_template.keys())
    shape_line="\t".join(str(1)*len(color_template))

    anno_out.write("LEGEND_COLORS\t"+color_line+"\n")
    anno_out.write("LEGEND_LABELS\t"+label_line+"\n")
    anno_out.write("LEGEND_SHAPES\t"+shape_line+"\n")
    anno_out.write("LEGEND_TITLE\tClade\n")
    anno_out.write("MARGIN\t5\n")
    anno_out.write("STRIP_WIDTH\t25\n")
    anno_out.write("DATA\n")

    tree=Tree(treefile)
    # gene_name=[]
    for i in tree:
        sub=define_sub(i.name,subg_clade_dic)
        c_line="\t".join([i.name,color_template[sub],sub])+"\n"
        anno_out.write(c_line)
    anno_out.close()
    print("Annotation successfully stored in "+treefile+"_itol_annotation.txt")

def main(treefile, config):
    print(f"Tree file: {treefile}")
    print(f"Config file: {config}")
    clade_subg_dic={}
    subg_clade_dic={}
    with open(config,'r') as c_file:
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

    cur_color = generate_random_colors(len(clade_subg_dic))
    color_template = dict(zip(clade_subg_dic.keys(), cur_color))

    write_itol_anno(treefile,template_color,subg_clade_dic)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate itol annotation file for AGFD tree.")
    parser.add_argument("-t","--treefile", type=str, help="Newick tree file used in AGFD.")
    parser.add_argument("-c","--config", type=str, help="Config file used in AGFD.")
    args = parser.parse_args()
    main(args.treefile, args.config)
    # test_treefile=r'E:\Bio_analysis\HGT_newcollection\HGTfinder_script\AGFD_V2_0910\Example_file\Testfile2.treefile'
    # config_file=r'E:\Bio_analysis\HGT_newcollection\HGTfinder_script\AGFD_V2_0910\Example_file\test_config.txt'
    # main(test_treefile,config_file)