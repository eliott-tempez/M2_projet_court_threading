## Double dynamic programming threading program
#### Eliott TEMPEZ - M2BI - Université Paris Cité

*Description*

All fasta and pdb files in the data/proteins/ folder were found on the PDB databate, and are classified into 3 fold classes : all alpha, all beta, and alpha and beta (a+b). They are all between 40 and 60 residues-long. There are 4 proteins per class, meaning the programs were tested on 12 different proteins.The data/small_prot/ folder contains a fasta and pdb file for a 10-residues long protein, which was used to develop the program.

To read the full report on this project, you can find a pdf file in the doc/ folder.


### How-to
#### Prerequisites (from the current directory, in a linux terminal):
* [install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
* install the conda environment
`conda env create -f environment.yaml`
* activate the conda environment
`conda activate projet_court_env`

#### To align a single protein sequence of interest to a template structure:
run the code line:
`python src/main.py arg1 arg2 (arg3)`
* arg1 : path to the fasta file of the sequence of interest
* arg2 : path to the pdb file of the template protein
* arg3 (optional) : gap penalty

#### To run the program on several templates for one sequence of interest:

