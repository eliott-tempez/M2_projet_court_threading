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
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import seq1


DOPE_FILE = "data/dope.par"
GAP_PENALTY = 0
#### A SUPPRIMER POUR AJOUTER LES COMMANDES UTILISATEURS
seq_file = "data/5AWL.fasta"
template_file = "data/5AWL.pdb"



##############################################################################
###################                CLASSES                ####################
##############################################################################

class AlphaCarbon:
    """Information on alpha carbons in a 3D structure"""
    def __init__(self, res_number, x, y, z):
        self.atom_name = str(res_number)
        # coordinates
        self.x = x
        self.y = y
        self.z = z
        
    def calculate_distance(self, carbon2):
        """Calculate the distance in angstroms between 2 alpha carbons"""
        return math.sqrt((carbon2.x-self.x)**2 + 
                         (carbon2.y-self.y)**2 + 
                         (carbon2.z-self.z)**2)
        
        
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
    
    def get_length(self):
        """Get number of alpha carbons in the protein."""
        return len(self.carbon_list)
    
    def calculate_inter_ca_distance(self, n_ca1, n_ca2):
        """Return distance between two alpha carbons in angstroms.

        Args:
            n_ca1 (int): number of the first carbon in order in the protein
            n_ca2 (int): number of the second carbon in order in the protein

        Returns:
            float: distance in angstroms
        """
        return self.carbon_list[n_ca1].calculate_distance(self.carbon_list[n_ca2])
        
    

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
    

def check_positive_number(number):
    """Check if a number is greater than zero"""
    assert number >= 0, f"Error : a distance cannot be negative, \
        check your pdb file"
    
    
def check_matrix_init_value(value):
    value_is_number = isinstance(value, (int, float))
    value_is_inf = np.isinf(value)
    assert value_is_number or value_is_inf, f"Error : matrix should \
        be initiated with a number or infinite value"



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
    # change 3-letter amino-acid code to 1-letter
    dope_mat_ca.loc[:, "res1"] = [seq1(res) for res in dope_mat_ca["res1"]]
    dope_mat_ca.loc[:, "res2"] = [seq1(res) for res in dope_mat_ca["res2"]]
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
    return seq_str


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


####----------------------      HANDLE MATRICES      ---------------------####

def get_dope_score(dope_mat, res1, res2, distance):
    """Get corresponding dope score from dope dataframe.

    Args:
        dope_mat (pd.DataFrame): the dope scores dataframe
        res1 (str): name of the first residue (3-letter in all caps)
        res2 (str): name of the second residue
        distance (float): distance between the two residues in angstroms

    Returns:
        float: corresponding dope score
    """
    check_positive_number(distance)
    # return high value if distance too short
    if distance < 0.25:
        return 10
    # return null value if distance too high
    elif distance > 14.75:
        return 0
    # if we have the exact distance in the dope dataframe, return value
    elif distance in dope_mat.columns:
        line_dope_mat = dope_mat[(dope_mat["res1"] == min(res1, res2)) & 
                        (dope_mat["res2"] == max(res1, res2))]
        return line_dope_mat[distance].iloc[0]
    # if we don't, apply an affine function
    else:
        line_dope_mat = dope_mat[(dope_mat["res1"] == min(res1, res2)) & 
                        (dope_mat["res2"] == max(res1, res2))]
        # extract distances in dataframe surrounding our distance value
        numeric_columns = [col for col in dope_mat.columns 
                           if isinstance(col, (int, float))]
        distance_below = max([x for x in numeric_columns if x < distance])
        distance_above = min([x for x in numeric_columns if x > distance])
        score_below = line_dope_mat[distance_below].iloc[0]
        score_above = line_dope_mat[distance_above].iloc[0]
        # return proportional dope score
        score_interpolated = score_below + (((distance - distance_below) / 
                                            (distance_above - distance_below)) * 
                                            (score_above - score_below))
        return round(score_interpolated, 2)


def initialise_matrix(shape, value=0.0):
    """Initialise numpy matrix of given shape with chosen value"""
    check_matrix_init_value(value)
    return np.full(shape, value)


def fill_distance_matrix(template_prot):
    """Create a distance matrix for all pairwise distances in the template."""
    n = template_prot.get_length()
    distance_matrix = initialise_matrix((n, n))
    # calculate all distances between alpha carbons in the protein
    for i in range(n):
        for j in range(n):
            if i !=j:
                dist = template_prot.calculate_inter_ca_distance(i, j)
                distance_matrix[i, j] = round(dist, 2)
    return distance_matrix


def fill_low_level_matrix(shape, dist_matrix, dope_matrix, 
                          test_sequence, gap_penalty, i, j):
    """Fill and return a singular low-level matrix for the start point [i, j].

    Args:
        shape (tuple): shape of the matrix
        dist_matrix (np.ndarray): distance matrix for the template
        dope_matrix (pd.DataFrame): dataframe of the dope scores
        test_sequence (list): list of residues in the test sequence
        i (int): line number of the start point for this matrix
        j (int): column number of the start point for this matrix
        gap_penalty (int): gap penalty

    Returns:
        np.ndarray: low-level matrix for the start point [i, j]
    """
    # fill matrix with inf because we won't fill some areas 
    # and the goal is to minimise
    L_mat = initialise_matrix(shape)
    n_query = shape[0]
    n_template = shape[1]
    
    L_mat[i, j] = 0
    # fill matrix area before the starting point 
    # but avoid the two residues next to starting point
    for k in range(i-3, -1, -1):
        for l in range(j-3, -1, -1):
            dist = dist_matrix[j, l]
            res1 = test_sequence[i]
            res2 = test_sequence[k]
            score = get_dope_score(dope_matrix, res1, res2, dist)
            # calculate minimum score
            L_mat[k, l] = min(
                L_mat[k+1, l+1] + score,
                L_mat[k, l+1] + gap_penalty,
                L_mat[k+1, l] + gap_penalty)
            
    # fill matrix area after the starting point
    for k in range(i+3, shape[0]):
        for l in range(j+3, shape[1]):
            dist = dist_matrix[j, l]
            res1 = test_sequence[i]
            res2 = test_sequence[k]
            score = get_dope_score(dope_matrix, res1, res2, dist)
            # calculate minimum score
            L_mat[k, l] = min(
                L_mat[k-1, l-1] + score,
                L_mat[k, l-1] + gap_penalty,
                L_mat[k-1, l] + gap_penalty)
    #return matrix
    return L_mat


def create_global_low_level(shape, dist_matrix, dope_matrix, 
                            test_sequence, gap_penalty):
    """Create the super-low-level matrix filled with each low-level matrix"""
    # initialise global matrix
    global_L_mat = np.empty(shape, dtype=object)
    # for each cell, create a low-level matrix
    for i in range(shape[0]):
        for j in range(shape[1]):
            L_mat = fill_low_level_matrix(shape, dist_matrix, dope_matrix,
                            test_sequence, gap_penalty, i, j)
            global_L_mat[i, j] = L_mat
    return global_L_mat
            
            
        
    
                    
                
            
        
    












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
    
    ####-------------    SET UP MATRICES    -------------####
    # distance matrix
    dist_matrix = fill_distance_matrix(template_prot)
    # low and high level matrices
    row_names_query = test_seq
    n_atoms_template = template_prot.get_length()
    mat_shape = (len(row_names_query), n_atoms_template)
    
    # create low-level-matrix
    global_L_mat = create_global_low_level(mat_shape, dist_matrix,
                                           dope_scores, test_seq, GAP_PENALTY)
    
    
    
    
    
    
    
    
    
    
    
    
    low_lvl_mat = initialise_matrix(mat_shape)
    high_lvl_mat = initialise_matrix(mat_shape)
    
    
    
    
    
    
    
    
    
            

