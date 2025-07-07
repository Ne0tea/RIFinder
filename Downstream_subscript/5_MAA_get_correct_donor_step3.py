'''
Descripttion: 解析kaks结果文件，更新其中异常kaks值
Author: Ne0tea
version: 
Date: 2023-11-02 15:16:43
LastEditors: Ne0tea
LastEditTime: 2023-11-03 15:43:04
'''

import sys

def main(kf1,kf2):
    new_kaks={}
    outfile=open('gene_adjusted_kaks.txt','w')
    with open(kf2,'r') as kk2:
        for i in kk2:
            line=i.strip()
            gene=line.split(' ')[0]
            new_kaks[gene]=line
    with open(kf1,'r') as kk1:
        for i in kk1:
            line=i.strip()
            gene=line.split(' ')[0]
            if gene in new_kaks:
                outfile.write(new_kaks[gene]+'\n')
            else:
                outfile.write(line+'\n')
    outfile.close()


if __name__ == "__main__":
    kaks_file1=sys.argv[1]
    kaks_file2=sys.argv[2]
    # kf1 old file | kf2 new file
    main(kaks_file1,kaks_file2)