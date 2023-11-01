# Horizontal_gene_transfer_project
Identify Horizontal gene transfer gene across Poaceae subfamily
usage: HgtFinder [-h] -t TREE -c CONFIG -go GO -o OUT [--kmax KMAX] [--kmin KMIN] [--support SUPPORT]
-------------------------------------------------------------------------------------------------------
HgtFinder
Author: Yujie Huang <12116008@zju.edu.cn>, ZJU
Version: v1.0
Identify HGT events based on gene tree
-------------------------------------------------------------------------------------------------------
options:
  -h, --help            show this help message and exit
  -t TREE, --tree TREE  supply a NWK gene tree file set
  -c CONFIG, --config CONFIG
                        supply config file indicating phlogenic relationship and clade name
  -go GO, --gene-order GO
                        supply gene order file
  -o OUT, --out OUT     supply out info file
  --kmax KMAX           the porpotion of dominate family in KNN cluster
  --kmin KMIN           the porpotion of minor family in KNN cluster
  --support SUPPORT     Filter out branches with bootstrap values lower than the standard during the identification of ancient transfer events.
