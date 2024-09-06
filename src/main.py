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
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO



DOPE_FILE = "data/dope.par"
#### A SUPPRIMER POUR AJOUTER LES COMMANDES UTILISATEURS
seq_file = "data/3NIR.fasta"
template_file = "data/3NIR.pdb"



#################################################################################
#####################                CLASSES                #####################
#################################################################################

class AlphaCarbon:
    """Information on alpha carbons in a 3D structure"""
    def __init__(self, residue_name, residue_number, x, y, z):
        self.residue_name = residue_name
        self.residue_number = residue_number
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
            residue (AphaCarbon): alpha carbon to add to the end of the existing protein
        """
        self.carbon_list.append(residue)
        
    def get_residues_list(self):
        """Get list of residues in the protein

        Returns:
            list: list of residues in 3-letters form
        """
        return [ca.residue_name for ca in self.carbon_list]
        
    







#################################################################################
####################                FUNCTIONS                ####################
#################################################################################

####------------------------      HANDLE ERRORS      ------------------------####
def check_file_exists(file_path):
    """Check if file exists"""
    assert os.path.isfile(file_path), f"Error : The file '{file_path}' doesn't exist"
    

def check_pdb_file(file_path):
    """Check if pdb file is conform"""
    assert file_path.split(".")[-1] == "pdb", f"Error : The file '{file_path}' is not a pdb file"
    
def check_protein(file_path, protein):
    """Check if protein contains atoms ; if not, the problem is the pdb file"""
    assert protein.carbon_list != [], f"Error : The pdb file '{file_path}' is empty or non valid"




####--------------------------      READ FILES      -------------------------####
def get_dtf_from_dope_file(file_path):
    """Read dope file and return pandas dataframe.

    Args:
        file_path (str): dope file name

    Returns:
        pd.DataFrame: dataframe of dope score for each pairwise interaction
    """
    # check if file existe
    check_file_exists(file_path)
    # get the names of the distances corresponding to the dope scores
    anstrom_len_as_str = np.arange(0.25, 15, 0.5).tolist()
    # read dope file
    dope_mat_full = pd.read_csv(file_path, sep=" ",
                                names=["res1", "atom1", "res2", "atom2"] + anstrom_len_as_str)
    # keep only carbon alpha dope scores
    dope_mat_ca = dope_mat_full[(dope_mat_full["atom1"] == "CA") & (dope_mat_full["atom2"] == "CA")]
    # return dataframe without unimportant columns
    return dope_mat_ca.drop(columns=["atom1", "atom2"]).reset_index(drop=True)


def get_query_from_file(file_path):
    """Read simple fasta file and return sequence.

    Args:
        filename (str): fasta file name

    Returns:
        str: query sequence
    """
    # check if file exists
    check_file_exists(file_path)
    # extract sequence
    return str(SeqIO.read(file_path, "fasta").seq)


def get_template_from_file(file_path):
    # check if file exists and is pdb file
    check_file_exists(file_path)
    check_pdb_file(file_path)
    # create protein
    template = Protein()
    # parse pdb file
    with open(file_path, "r") as pdb_file:
        for line in pdb_file:
            # check if line is an alpha carbon
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                # extract atom informations
                res_name = line[17:20].strip()
                res_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                alpha_c = AlphaCarbon(res_name, res_number, x, y, z)
                # add all alpha carbons to the template protein
                template.add_residue(alpha_c)
    # check protein
    check_protein(file_path, template)
    return template




if __name__ == "__main__":
    # read dope values
    dope_scores = get_dtf_from_dope_file(DOPE_FILE)
    # read query sequence
    test_seq = get_query_from_file(seq_file)
    # read 3D data
    template = get_template_from_file(template_file)
    # get template protein sequence
    template_sequence = template.get_residues_list()
    print(template_sequence)
    
            

