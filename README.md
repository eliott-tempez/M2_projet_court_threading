## Double dynamic programming threading program
#### Eliott TEMPEZ - M2BI - Université Paris Cité

This project aims to implement a threading program, a method used to model proteins which have the same fold as proteins of known structures, with python. It is based on "double" dynamic programming, as explained in Jones et al. (1995). You can find the scripts in the src/ directory, and the results in the results/ directory.

All fasta and pdb files in the data/proteins/ folder were found on the PDB databate, and are classified into 3 fold classes : all alpha, all beta, and alpha and beta (a+b). They are all between 40 and 50 residues-long. There are 3 proteins per class, meaning the programs were tested on 9 different proteins.The data/small_prot/ folder contains a fasta and pdb file for a 10-residues long protein, which was used to develop the program.

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
`python src/align_structure.py arg1 arg2 (arg3)`
* arg1 : path to the fasta file of the sequence of interest
* arg2 : path to the pdb file of the template protein
* arg3 (optional) : gap penalty

#### To run the program on several templates for one sequence of interest:
`python src/meta_alignment.py arg1 arg2`
* arg1 : directory path containing pdb files that can be in subdirectories or not
* arg2 : path to the fasta file of the protein of interest
This program returns a text file in results/alignments/ with all results from the main program *align_structure.py*

#### To run the program evaluating the alignment of a protein of interest on several templates:
`python src/evaluate_alignment.py`
This program automatically processes the files in the results/alignments/ directory and outputs graphs in the results/eval/ directory.
