'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-01-01 16:16:43
LastEditors: Ne0tea
LastEditTime: 2025-01-01 16:22:14
'''
from collections import defaultdict
import sys
def concatenate_fasta_sequences(fasta_file):
    # 字典用于存储以相同前缀的序列
    sequences = defaultdict(str)
    sp_seq_dic={}
    with open(fasta_file, 'r') as f:
        header = None
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):  # 序列的标题行
                if header is not None:
                    sequences[header.split('_')[0]] += sequence
                header = line[1:]  # 去掉 '>' 字符
                sequence = ""  # 重置序列
                if header.split('_')[0] not in sp_seq_dic:
                    sp_seq_dic[header.split('_')[0]]=1
                else:
                    sp_seq_dic[header.split('_')[0]]+=1
            else:
                sequence += line  # 将序列行拼接到一起

        # 处理最后一个序列
        if header is not None:
            sequences[header.split('_')[0]] += sequence
    print(sp_seq_dic)
    # for prefix, seq in sequences.items():
    #     print(f">{prefix}")
    #     print(seq)

# 调用函数，替换文件路径为你的FASTA文件路径
concatenate_fasta_sequences(sys.argv[1])
