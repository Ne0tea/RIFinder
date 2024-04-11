# Horizontal_gene_transfer_project  
Identify Horizontal gene transfer gene across Poaceae subfamily  

HgtFinder  
-------------------------------------------------------------------------------------------------------  

##</a>Getting Started
```sh
# Install HGTfinder which  is just a straightforward Python workflow(requiring sklearn and ete3)
git clone https://github.com/Ne0tea/HGTfinder

# Run on test data (use -f0 for small datasets)
python HGTfinder.py -t ./example_data/test_tree_set_file.txt -c ./example_data/test_config.txt -go ./example_data/test_gene_order.txt -o ./example_data/test_out_file.txt

# Run based on a set of gene tree
python HGTfinder.py -t TREE_SET -c CONFIG -go GO -o OUT

# help documentation
usage: HgtFinder [-h] -t TREE -c CONFIG -go GO -o OUT [--kmax KMAX] [--kmin KMIN] [--support SUPPORT]  
-------------------------------------------------------------------------------------------------------
options:  
  -h, --help            show this help message and exit  
  -t TREE, --tree TREE  supply a NWK gene tree file set  
  -c CONFIG, --config CONFIG  
                        supply config file indicating phlogenic relationship and clade name  
  -go GO, --gene-order GO  
                        supply sorted gene order file  
  -o OUT, --out OUT     supply out info file  
  --kmax KMAX           the porpotion of dominate family in KNN cluster  
  --kmin KMIN           the porpotion of minor family in KNN cluster  
  --support SUPPORT     Filter out branches with bootstrap values lower than the standard during the identification of ancient transfer events.  
```

## </a> issue
If any questions, contact at Yujie Huang <yujiehuang@zju.edu.cn>, or Dongya Wu <wudongya@zju.edu.cn>,ZJU  

## </a>Citating HgtFinder

If you use HgtFinder in your work, please wait:
......
