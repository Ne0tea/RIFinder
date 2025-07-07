'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2023-11-05 11:02:33
LastEditors: Ne0tea
LastEditTime: 2023-11-05 23:08:34
'''
from ete3 import Tree
import sys
import os

anno_file=sys.argv[1]


with open(anno_file,'r') as af:
    for i in af:
        line=i.strip()
        gene=line.split()[0]
        chr=line.split()[1]
        start=int(line.split()[4])
        end=int(line.split()[5])
        file_stat=os.system('ls /public/home/huangyj/grass_horizatal/ref/genome/'+str(gene[:4])+'.fasta')
        if file_stat:
            print(line[0],'there is no according fasta')
        else:
            os.system('seqkit subseq --chr '+chr+' -r '+str(start)+':'+str(end)+' -o '+gene+'.seq'\
                +' /public/home/huangyj/grass_horizatal/ref/genome/'+str(gene[:4])+'.fasta')
# Aspl_GWHBHQB00000023_3446002-3448235 Aspl_GWHBHQB00000023
            ori_name=chr+'_'+str(start)+'-'+str(end)+' '+chr
            os.system('sed -i \'s/'+ori_name+'/'+gene+'/g\' '+gene+'.seq')
print('Success!')