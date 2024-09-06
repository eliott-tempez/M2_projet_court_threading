"""This program is a threading program, a method used to model proteins 
which have the same fold as proteins of known structures. It is based
on "double" dynamic programming, as explained in Jones et al. (1995).

Usage:
======
    python main.py argument1 argument2

    argument1: a fasta file of the protein sequence you want to model
    argument2: a pdb file of the template protein
"""

__author__ = "Eliott TEMPEZ"
__date__ = "2024-09-11"


import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import seq3


DOPE_FILE = "data/dope.par"
#### A SUPPRIMER POUR AJOUTER LES COMMANDES UTILISATEURS
seq_file = "data/5AWL.fasta"
template_file = "data/5AWL.pdb"



##############################################################################
###################                CLASSES                ####################
##############################################################################

class AlphaCarbon:
    """Information on alpha carbons in a 3D structure"""
    def __init__(self, res_number, x, y, z):
        self.atom_name = "CA" + str(res_number)
        # coordinates
        self.x = x
        self.y = y
        self.z = z
        
        
class Protein:
    """Protein class which includes all alpha carbons."""
    def __init__(self):
        # list of alpha carbons 
        self.carbon_list = []
    
    def add_residue(self, residue):
        """Add a residue to the protein.

        Args:
            residue (AphaCarbon): carbon to add at the end of the protein
        """
        self.carbon_list.append(residue)
    
    def get_carbons_list(self):
        """Get list of atoms (one per residue) in the protein

        Returns:
            list: list of atom names
        """
        return [ca.atom_name for ca in self.carbon_list]
        
    

##############################################################################
###################                FUNCTIONS                ##################
##############################################################################

####-----------------------      HANDLE ERRORS      ----------------------####
def check_file_exists(file_path):
    """Check if file exists"""
    fils_exists = os.path.isfile(file_path)
    assert fils_exists, f"Error : The file '{file_path}' doesn't exist"
    

def check_pdb_file(file_path):
    """Check if pdb file is conform"""
    is_pdb = file_path.split(".")[-1] == "pdb"
    assert is_pdb, f"Error : The file '{file_path}' is not a pdb file"
    

def check_fasta_file(file_path):
    """Check if fasta file is conform"""
    is_fasta = file_path.split(".")[-1] in ["fasta", "fa"]
    assert is_fasta, f"Error : The file '{file_path}' is not a fasta file"


def check_sequence(file_path, seq):
    """Check if sequence contains residues:
       if not, the problem is the fasta file"""
    seq_not_empty = seq != ""
    assert seq_not_empty, f"Error : The fasta file '{file_path}' is empty"


def check_protein(file_path, protein):
    """Check if protein contains atoms: 
       if not, the problem is the pdb file"""
    prot_not_empty = protein.carbon_list != []
    assert prot_not_empty, f"Error : The pdb file '{file_path}' is not valid"


####------------------------      READ FILES      ------------------------####
def get_dtf_from_dope_file(file_path):
    """Read dope file and return pandas dataframe.

    Args:
        file_path (str): dope file path

    Returns:
        pd.DataFrame: dataframe of dope score for each pairwise interaction
    """
    # check if file existe
    check_file_exists(file_path)
    # get the names of the distances corresponding to the dope scores
    anstrom_len_as_str = np.arange(0.25, 15, 0.5).tolist()
    # read dope file
    colnames = ["res1", "atom1", "res2", "atom2"] + anstrom_len_as_str
    dope_mat_full = pd.read_csv(file_path, sep=" ", names=colnames)
    # keep only carbon alpha dope scores
    mask = (dope_mat_full["atom1"] == "CA") & (dope_mat_full["atom2"] == "CA")
    dope_mat_ca = dope_mat_full[mask]
    # return dataframe without unimportant columns
    return dope_mat_ca.drop(columns=["atom1", "atom2"]).reset_index(drop=True)


def get_query_from_fasta(file_path):
    """Read simple fasta file and return sequence.

    Args:
        filename (str): fasta file path

    Returns:
        str: query sequence
    """
    # check if file exists and is fasta file
    check_file_exists(file_path)
    check_fasta_file(file_path)
    # extract sequence in str form
    seq_str = str(SeqIO.read(file_path, "fasta").seq)
    check_sequence(file_path, seq_str)
    # transform sequence in 3-letter list form
    return [seq3(letter).upper() for letter in seq_str]


def get_template_from_pdb(file_path):
    """Read pdb file and extract information on all alpha carbons present.

    Args:
        file_path (str): pdb file path

    Returns:
        Protein: a Protein instance which contains infos on all alpha carbons
    """
    # check if file exists and is pdb file
    check_file_exists(file_path)
    check_pdb_file(file_path)
    # create empty protein
    template = Protein()
    # parse pdb file
    res_nb = 0
    with open(file_path, "r") as pdb_file:
        for line in pdb_file:
            # check if line is an alpha carbon
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                # check if we don't already have a position for this residue
                if res_nb != int(line[22:26].strip()):
                    res_nb += 1
                    # extract atom information
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    alpha_c = AlphaCarbon(res_nb, x, y, z)
                    # add all alpha carbons to the template protein
                    template.add_residue(alpha_c)
    # check protein
    check_protein(file_path, template)
    return template


####----------------------      HANDLE MATRIXES      ---------------------####

def initialise_matrix(shape):
    """Initialise numpy matrix of given shape with zeros"""
    return np.zeros(shape)












##############################################################################
#####################                MAIN                #####################
##############################################################################

if __name__ == "__main__":
    ####-------------       GET DATA       -------------####
    # read dope values
    dope_scores = get_dtf_from_dope_file(DOPE_FILE)
    # read query sequence
    test_seq = get_query_from_fasta(seq_file)
    # read 3D data
    template_prot = get_template_from_pdb(template_file)
    
    ####-------------    BUILD MATRIXES    -------------####
    row_names_query = test_seq
    col_names_template = template_prot.get_carbons_list()
    mat_shape = (len(row_names_query), len(col_names_template))
    low_lvl_mat = initialise_matrix(mat_shape)
    high_lvl_mat = initialise_matrix(mat_shape)
    
            

