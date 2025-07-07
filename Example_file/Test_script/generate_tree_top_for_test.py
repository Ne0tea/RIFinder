'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-14 15:55:39
LastEditors: Ne0tea
LastEditTime: 2024-09-18 14:06:16
'''
from itertools import combinations
from Inconsistant_toplogy_identify import decide_toplogy_from_process
import copy
import random
from ete3 import Tree
def standardize_tree_structure(tree):
    """
    对树进行标准化，使左右子树按照字典序排列，从而去除等价的树结构。
    """
    if tree.is_leaf():
        return tree
    
    # 对子树进行排序，如果有两个子节点
    if len(tree.children) == 2:
        left_child = standardize_tree_structure(tree.children[0])
        right_child = standardize_tree_structure(tree.children[1])
        
        # 获取左右子树的叶节点名称并进行排序比较
        left_names = sorted([leaf.name for leaf in left_child.iter_leaves()])
        right_names = sorted([leaf.name for leaf in right_child.iter_leaves()])
        
        if left_names > right_names:
            # 如果左子树的名称排序比右子树大，交换它们
            tree.swap_children()
    
    return tree

def generate_strict_binary_trees(nodes):
    if len(nodes) == 0:
        return []
    
    if len(nodes) == 1:
        tree = Tree()
        tree.name = nodes[0]
        return [tree]
    
    all_trees = []
    
    # 遍历每一个可能的左右子树组合（非相邻划分）
    for i in range(1, len(nodes)):
        # 使用 combinations 从 nodes 中挑选 i 个节点作为左子树
        for left_comb in combinations(nodes, i):
            left_comb = list(left_comb)
            right_comb = [n for n in nodes if n not in left_comb]
            
            # 递归生成所有可能的左子树和右子树
            left_subtrees = generate_strict_binary_trees(left_comb)
            right_subtrees = generate_strict_binary_trees(right_comb)
            
            # 合并左右子树，创建新的二叉树
            for left_subtree in left_subtrees:
                for right_subtree in right_subtrees:
                    root = Tree()
                    root.name = 'ROOT'
                    root.add_child(left_subtree)
                    root.add_child(right_subtree)
                    
                    # 标准化树结构
                    standardized_tree = standardize_tree_structure(root.copy())
                    all_trees.append(standardized_tree)

    return all_trees

def remove_duplicate_trees(trees):
    """
    去除等价的二叉树，通过将树的结构转换为字符串表示，去重。
    """
    unique_trees = []
    seen_structures = set()
    
    for tree in trees:
        # 将树结构转换为标准化字符串
        tree_str = tree.write(format=9)  # format=9 生成 Newick 树结构，不包含节点长度
        if tree_str not in seen_structures:
            seen_structures.add(tree_str)
            unique_trees.append(tree)
    
    return unique_trees

def get_sister_component(node):
    cur_leafs=[]
    cur_leafs.extend(node.get_leaf_names())
    cur_leafs.extend(node.get_sisters()[0].get_leaf_names())
    # print(cur_leafs)
    return cur_leafs
# 示例节点列表
subfamily_list=[['CHL', 'PAN'],['POO', 'ORY', 'BAM']]

clade_template_dic={'BAM': 7, 'ORY': 3, 'POO': 7, 'PAN': 8, 'CHL': 2, 'OUT': 1}
nodes = ['CHL', 'PAN','POO','ORY']
# 生成所有可能的严格二叉树
all_trees = generate_strict_binary_trees(nodes)
# 去除等价的树结构
unique_trees = remove_duplicate_trees(all_trees)
print(len(unique_trees))

# for tree in random.sample(unique_trees, 30):
# for tree in unique_trees:
#     for node in tree:
#         cur_process=[]
#         cur_node=copy.deepcopy(node)
#         cur_process.append([node.name])
#         # print(node.up,tree.write())
#         while cur_node.up:
#             cur_process.append(get_sister_component(cur_node))
#             cur_node=cur_node.up
#         a,b,c=decide_toplogy_from_process(cur_process,subfamily_list,'OUT',clade_template_dic)
#         if a == 'AGF':
#             print(cur_process)
#             print(tree)
#             print(a,b,c)
