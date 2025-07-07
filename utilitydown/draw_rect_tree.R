##########################################################################
# Script Name: draw_rect_tree.R
# 
# Description: 
# This script is designed to visualize a gene tree using the `ggtree` package. 
# It reads a phylogenetic tree from a specified file, roots the tree using 
# an outgroup, and colors the tree tips based on predefined species groups.
# The final tree is saved as a PDF file.
#
# The script uses various R libraries including ggtree, randomcoloR, 
# patchwork, and tidytree for tree manipulation and visualization.
#
# INPUT:
#  - tree_file: Path to the input Newick format tree file containing the gene tree.
#  - config_file: Path to the config file containing the phylogeny information, same as used in AGFD.
#  - color_template: A named list specifying the colors for each species group.
#  - tree_out_name: Path where the final rooted and annotated tree will be saved 
#    as a PDF.
# OUTPUT:
#  - tree_out_name
#
# Usage:
#  - This script can be run in an R environment by sourcing the script or running 
#    line by line. Ensure that all required packages are installed and the input 
#    tree file is available at the specified location.
#
# Example:
#  - Rscript draw_rect_tree.R <tree_file> <config_file> <tree_out_name>
#
# Libraries used:
#  - ggplot2
#  - ggtree
#  - randomcoloR
#  - treeio
#  - patchwork
#  - tidytree
#
# Author: Yujie Huang
# Date: 2024.10.17
##########################################################################

###INPUT need modified
# tree_file <- "E:/Bio_analysis/HGT_newcollection/4_example_recheck/GF_tree_set/GF_1881_intro_concatenate.contree"
# tree_out_name <- 'E:/Bio_analysis/HGT_newcollection/4_example_recheck/GF_tree_set/GF_1881_intro_concatenate.contree.pdf'
# config_file <- 'E:/Bio_analysis/HGT_newcollection/2_HGTfinder_workpath/Poaceae_config.txt'
###--------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
tree_file <- as.character(args[1])
config_file <- as.character(args[2])
tree_out_name <- as.character(args[3])
###--------------------------------------------------------------------------------------------

print("INPUT is plotting by ggtree.")
if (TRUE) {
  config_content <- readLines(config_file)
  subg_start <- grep(">subg list", config_content) + 1
  subg_lines <- config_content[subg_start:length(config_content)]
  subg_dic <- list()
  for (line in subg_lines) {
    # 以制表符分割
    parts <- strsplit(line, " ")[[1]]
    if (length(parts) > 1) {
      # 将第一个部分作为名称，其余部分作为值
      name <- parts[1]
      values <- parts[-1]
      subg_dic[[name]] <- values
    }
  }
  color_template<-c("BAM"="#b71515","ORY"="#e97a0c","POO"="#ffde0a","PAN"="#023e7d","CHL"="#a3cef1","OUT"="#a9a29c")
}

define_sub<-function(gene){
  sp=strsplit(gene,"_")[[1]][1]
  clade<-names(subg_dic)[sapply(subg_dic, function(x) sp %in% x)]
  return(clade)
}

suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(patchwork))
suppressMessages(library(tidytree))
tree<-read.tree(tree_file)

set_out=FALSE
sub_list=list(BAM=c(),ORY=c(),POO=c(),PAN=c(),CHL=c(),OUT=c())
gene_col<-c()
for (i in tree$tip.label) {
  sub<-define_sub(i)
  sub_list[[sub]]<-append(sub_list[[sub]],i)
  gene_col<-append(gene_col,color_template[sub])
  if (sub == "OUT" && set_out==FALSE) {
    tree <-root(tree, outgroup = i, edgelabel = TRUE)
    set_out=TRUE
  }
}
if (! set_out) {
  print("There is no OUT sp!")
}

tree_grouped<-groupOTU(tree,sub_list,"Sub_family")
tregraph_data<-ggtree(tree_grouped, layout="rectangular",  size=0.8)
tregraph<-tregraph_data+
  geom_tiplab(hjust=-0.05,size=4, color="white",geom = 'label', 
              label.size=0,fill=gene_col,fontface=2) +
  geom_tippoint(size=2)+
  geom_nodelab(geom="shadowtext")+
  theme(legend.position='none')+ geom_rootedge(rootedge = 0.01)+
  # geom_nodelab(geom="shadowtext",size=2)
  # geom_nodepoint(color="orange", alpha=1/4, size=2) +
  # theme_tree2()
  xlim(NA, max(tregraph_data$data$x)*1.8)

p_hg<-ceiling(length(gene_col)/5.5)
ggsave(tree_out_name,tregraph,height = p_hg,width=12,limitsize = FALSE)

print("Plotted successfully.")
