# Ancient Gene Flow Detection

## Overview
AGFD is a tool designed to identify gene flow events across a phylogeny of multiple assemblies. It analyzes gene tree structures and interspecies relationships to detect potential gene flow.

## Installation and Dependencies
Ensure the following Python libraries are installed:
```
ete3==3.1.3
matplotlib==3.7.1
numpy==1.23.5
pandas==2.2.3
scipy==1.10.1
seaborn==0.13.2
statsmodels==0.14.3
```
**To install any missing dependencies, you can use pip:**
`pip install -r requirements.txt`


## Running the Tool

### Command Line Arguments:
```
python agfd.py -f <gene_tree_file> -c <config_file> [-o <output_file>] [-rn] [-t <num_threads>] [--kmax <float>] [--kmin <float>] [--support <int>]
**Parameters:**

-f, --treefile: Gene tree file set, one tree file path per line.
-c, --config: Configuration file indicating phylogenetic relationships and subgroup species.
-o, --out: Output file name (default: AGFD_output.txt).
-rn, --rootnested: Whether to disentangle tree structure nested with Outgroup.
-t, --thread: Number of threads (default: 1).
--support: Filter out branches with bootstrap values lower than the standard during the identification of ancient transfer events (default: 70).
```
### Example Command:

`python agfd.py -f trees.txt -c config.txt -o output.txt -t 4 --support 70`


### Usage Example
Prepare your input files (`trees.txt`, `config.txt`).

#### Input Format

1. **trees.txt**: File containing the path to the gene tree file used to detect gene flow.
    ```
    /your/path/to/tree/file1.treefile
    /your/path/to/tree/file2.treefile
    /your/path/to/tree/file3.treefile
    ```

2. **config.txt**: Configuration of homoeologous chromosomes.
    ```
    >phylogeny							### > indicate the start of nwk phylogeny line, outgroup should be named as [OUT, outgroup];
    ((A,B),((C,D),E));					### suitable phylogeny format for func(Tree) in ete3;
    >subg list							### > indicate the clade of which species belong to;
    C C1 C2 C3 C4 C5 C6 				### all species abbreviations should be used as the prefix of the gene name in the trees; 
    D D1 D2 D3 D4 D5 D6 D7 D8
    E E1 E2 E3 E4 E5 E6 E7 E8
    A A1 A2 A3 A4 A5
    B B1 B2 B3 B4 B5
    OUT O1 O2
    ```

#### Output Format

1. **$prefix[**OUT**].out**: Main results file for tree topology conflicts. The file does not contain header files, but is output in the same order as in the example.

| Column Header     | Description                                                                                                                                                                                  | Example Value                                 |
|:-----------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:---------------------------------------------:|
| `GF id`           | Gene flow events based on the id of the detection order, one id per event.                                                                                                               | GF_1                                          |
| `OG id`           | Name of the tree file where the gene flow event occurred.                                                                                                                                 | OG0000153                                     |
| `GF gene number`  | Number of transferred genes involved in gene flow events.                                                                                                                                | 10                                            |
| `Conflict type`   | Labelling different types of conflicts between gene tree structure and species tree structure.                                                                                            | Topology inconsistent, Topology consistent, RGF |
| `Donor clade`     | Donor clade to which donor species belong in the gene flow event. For ancient gene flow events, the donor clade may be a coalescent unit of multiple existing clades, split by **|**. | A or A|B                                      |
| `Receive clade`   | Clade to which receiving species belong in the gene flow event. For ancient gene flow, a way consistent with the `donor clade` column will be used.                                        | A or A|B                                      |
| `Donor species`   | Donor species in the gene flow event. **|** was used to connect different species abbreviations, specified in the config file.                                                          | D1, D2|D3|D4                                  |
| `Receive species` | Receiving species in the gene flow event. The same method was used to connect different species as `Donor species`.                                                                       | A1, A2|A3|A4                                  |

2. **config.txt**: Configuration of homoeologous chromosomes.
    ```
    >phylogeny							### > indicate the start of nwk phylogeny line, outgroup should be named as [OUT, outgroup];
    ((A,B),((C,D),E));					### suitable phylogeny format for func(Tree) in ete3;
    >subg list							### > indicate the clade of which species belong to;
    C C1 C2 C3 C4 C5 C6 				### all species abbreviations should be used as the prefix of the gene name in the trees; 
    D D1 D2 D3 D4 D5 D6 D7 D8
    E E1 E2 E3 E4 E5 E6 E7 E8
    A A1 A2 A3 A4 A5
    B B1 B2 B3 B4 B5
    OUT O1 O2
    ```

Run the script using the command provided above. Check the output files for results.

#### Contact
For any issues or questions, please contact:
- Yujie Huang, yujiehuang@zju.edu.cn
- Dongya Wu, wudongya@zju.edu.cn

This README provides a comprehensive guide on how to install, run, and use the AGFD tool.
