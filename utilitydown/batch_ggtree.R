##########################################################################
# Script Name: batch_ggtree.R
# 
# Description: 
# This script is used to visualize phylogenetic trees with annotated clades and 
# outgroup rooting using the `ggtree` package. The script processes an input tree 
# file, groups species into clades based on a configuration file, and colors the 
# tree tips according to the clade. The output includes two tree visualizations: 
# one with branch lengths and one without branch lengths.
#
# INPUT:
#  - tree_file: Path to the Newick format tree file containing the phylogenetic tree.
#  - config_file: Path to the config file containing the phylogeny information, same as used in AGFD.
#  - tree_out_prefix: Prefix of the final tree visualizations (PDF) will be saved.
#
# OUTPUT:
#  - Two PDF files:
#    1. A tree with branch lengths preserved, saved as "<tree_out_prefix>_branchlength.pdf".
#    2. A tree without branch lengths (collapsed branches), saved as "<tree_out_prefix>.pdf".
#
# Notes:
#  - The script defines clades using species groups provided in the config file.
#  - The script roots the tree using the first species found in the 'OUT' group as the outgroup.
#  - If no 'OUT' species is found, the tree will not be rooted.
#
# Example:
#  Rscript batch_ggtree.R <tree_file> <config_file> <tree_out_prefix>
#
# Libraries used:
#  - ggplot2
#  - ggtree
#  - treeio
#  - patchwork
#
# Author: Yujie Huang
# Date: 2024.10.17
##########################################################################

###INPUT need modified
args <- commandArgs(trailingOnly = TRUE)
tree_file<-args[1]
config_file<-args[2]
tree_out_prefix<-args[3]
color_template=c("BAM"="#b71515","ORY"="#e97a0c","POO"="#ffde0a","PAN"="#023e7d","CHL"="#a3cef1","OUT"="#a9a29c")
###--------------------------------------------------------------------------------------------

print("INPUT is plotting by ggtree.")
if (TRUE) {
  config_content <- readLines(config_file)
  subg_start <- grep(">subg list", config_content) + 1
  subg_lines <- config_content[subg_start:length(config_content)]
  subg_dic <- list()
  for (line in subg_lines) {
    parts <- strsplit(line, " ")[[1]]
    if (length(parts) > 1) {
      name <- parts[1]
      values <- parts[-1]
      subg_dic[[name]] <- values
    }
  }  
}

define_sub<-function(gene){
  sp=strsplit(gene,"_")[[1]][1]
  clade<-names(subg_dic)[sapply(subg_dic, function(x) sp %in% x)]
  return(clade)
}

options (warn = -1)
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(patchwork))
tree<-read.tree(tree_file)

set_out=FALSE
sub_list=list(BB=c(),OO=c(),PP=c(),PAC=c(),MAD=c(),OUT=c())
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
tree_file_length<-ggtree(tree_grouped, layout="rectangular", size=0.8,ladderize=FALSE)
tregraph_length<-tree_file_length+
  geom_tiplab(hjust=-0.05,size=4, color="white",geom = 'label', 
              label.size=0,fill=gene_col,fontface=2) +
  geom_tippoint(size=2)+
  geom_nodelab(geom="shadowtext")+
  # scaleClade(NULL,scale=.1) +
  # geom_taxalink(data=trans_dat,
  #               mapping=aes(taxa1=gene_f,taxa2=gene_t,color=trans_color),
  #               size=2,ncp=3,offset=0.15)+
  theme(legend.position='none')+ geom_rootedge(rootedge = 0.01)+
  xlim(NA, max(tree_file_length$data$x)*1.8)

tree_file_nolength<-ggtree(tree_grouped, layout="rectangular", branch.length="none",size=0.8,ladderize=FALSE)
tregraph_nolength<-tree_file_nolength+
  geom_tiplab(hjust=-0.05,size=4, color="white",geom = 'label', 
              label.size=0,fill=gene_col,fontface=2) +
  geom_tippoint(size=2)+
  geom_nodelab(geom="shadowtext")+
  # scaleClade(NULL,scale=.1) +
  # geom_taxalink(data=trans_dat,
  #               mapping=aes(taxa1=gene_f,taxa2=gene_t,color=trans_color),
  #               size=2,ncp=3,offset=0.15)+
  theme(legend.position='none')+ geom_rootedge(rootedge = 0.01)+
  xlim(NA, max(tree_file_nolength$data$x)*1.8)

p_hg<-ceiling(length(gene_col)/5.5)
ggsave(paste(tree_out_prefix,"_branchlength.pdf",sep = ""),tregraph_length,height = p_hg,limitsize = FALSE)
ggsave(paste(tree_out_prefix,".pdf",sep = ""),tregraph_nolength,height = p_hg,limitsize = FALSE)

print("plotted successfully.")