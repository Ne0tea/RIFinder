'''
Descripttion:
Author: Ne0tea
version:
Date: 2023-08-27 23:23:23
LastEditors: Ne0tea
LastEditTime: 2023-10-15 16:49:43
'''
import sys
import difflib
import re
import os

def get_exon_size(exon_file:str,pre="out") -> dict:
    exon_l={}
    with open(exon_file, 'r') as f:
        for line in f:
            lin = line.strip().split()
            a = len(lin)
            for i in range(6, a, 2):
                exon = abs(int(lin[i]) - int(lin[i-1])) + 1
                exon_l[lin[0]]=exon
    return exon_l
#上述程序则是计算每个外显子的大小

def get_exon_num(exon_file:str,pre="out") -> dict:
    exon_n={}
    with open(exon_file, 'r') as f:
        for line in f:
            lin = line.strip().split()
            a = len(lin)
            n = (a - 5)/2
            exon_n[lin[0]]=n
    return exon_n

#上述这个程序则是统计每个基因外显子的数量

def get_intro_size(exon_file:str,pre="out") -> dict:
    intro={}
    with open(exon_file, 'r') as f:
        for line in f:
            intro_size=0
            intro_num=0
            lin = line.strip().split()
            a = len(lin)
            if a == 7:
                intro[lin[0]]=0
            if a > 7:
                if lin[4] == '+':
                    for i in range(7, a, 2):
                        intron = abs(int(lin[i]) - int(lin[i-1]) - 1)
                        intro_size+=intron
                        intro_num+=1
                        # intro_l[lin[0]]=intron
                if lin[4] == '-':
                    for i in range(8, a, 2):
                        intron = abs(int(lin[i]) + 1 - int(lin[i - 3]))
                        intro_num+=1
                        intro_size+=intron
                        # intro_l[lin[0]]=intron
            intro[lin[0]]=[intro_num,intro_size]
    #intro['gene1']=[3,1000]
    return intro
#上述这个程序则是计算所有内含子的大小

def main(file:str,sca_len_file,intro=True,gsize=False,check=False,static_mode=False):
    p='(\w+)_full'
    sp=re.findall(p,file)[0]
    exon_file=open('tmp_exon.txt','w')
    type=[]
    sca_len={}

    with open(sca_len_file,'r') as len_f:
        for i in len_f:
            line=i.strip().split()
            sca_len[line[0]]=int(line[1])

    with open(file,'r') as pre_f:
        lines=[line for line in pre_f if not line.startswith('#')][:100]
        for line in lines:
            if not line.strip():
                continue
            lin = line.strip().split()
            type.append(lin[2])
        if type.count('CDS') == 0 and type.count('exon') != 0:
            count_type='exon'
        elif type.count('CDS') != 0 and type.count('exon') == 0:
            count_type='CDS'
        elif type.count('CDS') != 0 and type.count('exon') != 0:
            count_type='exon'
        else:
            print('Error! no suitable count type in ',os.path.basename(file))
            return
        if type.count('gene') == 0 and type.count('mRNA') != 0:
            name_type='mRNA'
        elif type.count('gene') != 0 and type.count('mRNA') == 0:
            name_type='gene'
        elif type.count('CDS') != 0 and type.count('mRNA') != 0:
            name_type='gene'
        else:
            print('Error! no suitable name type in ',os.path.basename(file))
            return
    with open(file, 'r') as f:
        gene_l={}
        cds_l=[]
        mRNA_l=[]
        num=0
        for line in f:
            # print(line)
            if line.startswith('#') or not line.strip():
                continue
            lin = line.strip().split()
            # print(lin[8].split(';')[1])
            if lin[2] == name_type:
                name = lin[8].split(';')[0].split('=')[-1]
                if num != 0:
                    exon_file.write("\n")
                exon_file.write(" ".join([name, lin[0], lin[3], lin[4],lin[6]])+' ')
                # print(line[0])
                # cur_chr=difflib.get_close_matches(str(sp)+str(line[0]),sca_len.keys(),2, cutoff=0.4)[0]
                # print(str(lin[0]))
                cur_chr=str(sp)+"_"+str(lin[0])
                if cur_chr not in sca_len:
                    cur_chr=difflib.get_close_matches(str(sp)+"_"+str(line[0]),sca_len.keys(),1, cutoff=0.4)[0]
                    gene_l[name]=[abs(int(lin[3])-int(lin[4])),round(int(lin[3])/sca_len[cur_chr],2)]
                else:
                    gene_l[name]=[abs(int(lin[3])-int(lin[4])),round(int(lin[3])/sca_len[cur_chr],2)]
                num+=1
                # print(name)
            if lin[2] == count_type:
                exon_file.write(" ".join([lin[3], lin[4]])+" ")
            if type.count('CDS') != 0 and lin[2] == "CDS":
                cds_l.append(abs(int(lin[3])-int(lin[4])))
            if type.count('mRNA') != 0 and lin[2] == "mRNA":
                mRNA_l.append(abs(int(lin[3])-int(lin[4])))
    exon_file.close()
    # cds_num = len(open('tmp_exon.txt','r').readlines())
    size_dic=get_exon_size('tmp_exon.txt')
    num_dic=get_exon_num('tmp_exon.txt')
    intro_dic=get_intro_size('tmp_exon.txt')

    if check:
        sizeF=next(iter(size_dic))
        sizeL=size_dic.get(next(iter(size_dic)))
        # print(sizeF,sizeL)
        numF=next(iter(num_dic))
        numL=num_dic.get(next(iter(num_dic)))
        intrF=next(iter(intro_dic))
        intrL=intro_dic.get(next(iter(intro_dic)))
        print(sizeF,sizeL,numF,numL,intrF,intrL)
    if static_mode:
        print('gene number is :',num)
        print('mRNA number is :',len(mRNA_l))
        print('CDS number is :',len(cds_l))
        print('max gene length is :',max(gene_l.values()))
        print('mean gene length is :',sum(gene_l.values())/len(gene_l))
        print('max CDS length is :',max(cds_l))
        print('mean CDS length is :',sum(cds_l)/len(cds_l))
        print('max exon length is :',max(size_dic.values()))
        print('mean exon length is :',sum(size_dic.values())/sum(num_dic.values()))
        print('mean exon num per CDS is :',sum(num_dic.values())/len(cds_l))
        print('max intro length is :',max(intro_dic.values()))
        print('mean intro length is :',sum(intro_dic.values())/sum(num_dic.values()))
    # print(intro_dic['Ola018545.1'])
    if intro and gsize:
        return intro_dic,gene_l
    elif intro and not gsize:
        return intro_dic




if __name__ == "__main__":
    # gff3_file=sys.argv[1]
    # sca_len_file=sys.argv[2]
    gff3_file=r"E:\Bio_analysis\Grass_horizatal_transversion\8_make_additional\prase_gff\Cson_full.gff3"
    sca_len_file=r"E:\Bio_analysis\Grass_horizatal_transversion\8_make_additional\genome.size"
    intro_dic,gene_dic=main(gff3_file,sca_len_file,gsize=True,check=False,static_mode=False)
    print(intro_dic,gene_dic)