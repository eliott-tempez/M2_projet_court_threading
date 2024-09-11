"""
This program is a threading program, a method used to model proteins
which have the same fold as proteins of known structures. It is based
on "double" dynamic programming, as explained in Jones et al. (1995).

Usage:
======
    python main.py argument1 argument2

    argument1: a fasta file of the protein sequence you want to model
    argument2: a pdb file of the template protein
    argument3 (optional) : the gap penalty for the alignment, default value 0

    Returns :
        The optimal alignment of the sequence and the template protein
        The total energy of this alignment
"""

__author__ = "Eliott TEMPEZ"
__date__ = "2024-09-11"


import os
import argparse
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import seq1


DOPE_FILE = "data/dope.par"


##############################################################################
###################                CLASSES                ####################
##############################################################################

class AlphaCarbon:
    """Information on alpha carbons in a 3D structure."""

    def __init__(self, res_number, x, y, z):
        """Alpha carbon constructor.

        Args:
            res_number (int): residue number (from 1)
            x (int): x coordinate
            y (int): y coordinate
            z (int): z coordinate
        """
        self.atom_name = str(res_number)
        # coordinates
        self.x = x
        self.y = y
        self.z = z

    def calculate_distance(self, carbon2):
        """Calculate the distance in angstroms between 2 alpha carbons."""
        return math.sqrt((carbon2.x-self.x)**2 +
                         (carbon2.y-self.y)**2 +
                         (carbon2.z-self.z)**2)


class Protein:
    """Protein class which includes all alpha carbons."""

    def __init__(self):
        """Protein constructor."""
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
        c1 = self.carbon_list[n_ca1]
        c2 = self.carbon_list[n_ca2]
        return c1.calculate_distance(c2)


##############################################################################
###################                FUNCTIONS                ##################
##############################################################################

####-----------------------      HANDLE ERRORS      ----------------------####
def check_file_exists(file_path):
    """Check if file exists."""
    fils_exists = os.path.isfile(file_path)
    assert fils_exists, f"Error : The file '{file_path}' doesn't exist"


def check_pdb_file(file_path):
    """Check if pdb file is conform."""
    is_pdb = file_path.split(".")[-1] == "pdb"
    assert is_pdb, f"Error : The file '{file_path}' is not a pdb file"


def check_fasta_file(file_path):
    """Check if fasta file is conform."""
    is_fasta = file_path.split(".")[-1] in ["fasta", "fa"]
    assert is_fasta, f"Error : The file '{file_path}' is not a fasta file"


def check_sequence(file_path, seq):
    """Check if sequence contains residues.
    If not, the problem is the fasta file
    """
    seq_not_empty = seq != ""
    assert seq_not_empty, f"Error : The fasta file '{file_path}' is empty"


def check_protein(file_path, protein):
    """Check if protein contains atoms.
    If not, the problem is the pdb file
    """
    prot_not_empty = protein.carbon_list != []
    assert prot_not_empty, f"Error : The pdb file '{file_path}' is not valid"


def check_positive_number(number):
    """Check if a number is greater than zero."""
    assert number >= 0, f"Error : a distance cannot be negative, \
        check your pdb file"


####------------------------      READ ARGS      -------------------------####
def read_args():
    """Read and return command line arguments.

    Returns:
        argparse.Namespace: command lign arguments
    """
    descr = "Threading program for protein modeling \
        based on double dynamic programming."
    parser = argparse.ArgumentParser(description=descr)

    # arguments
    fasta_descr = "fasta file of the protein sequence you want to model"
    pdb_descr = "pdb file of the template protein"
    gap_descr = "gap penalty for the alignment, default value 0"

    parser.add_argument("seq_file", type=str, help=fasta_descr)
    parser.add_argument("template_file", type=str, help=pdb_descr)
    parser.add_argument("gap_penalty", type=float,
                        nargs='?', default=0, help=gap_descr)
    args = parser.parse_args()
    return args


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
    anstrom_len = np.arange(0.25, 15, 0.5).tolist()
    # read dope file
    colnames = ["res1", "atom1", "res2", "atom2"] + anstrom_len
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
    # check sequence isn't empty
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
                    # create Atom instance
                    alpha_c = AlphaCarbon(res_nb, x, y, z)
                    # add all alpha carbons to the template protein
                    template.add_residue(alpha_c)
    # check protein isn't empty
    check_protein(file_path, template)
    return template


####----------------------      HANDLE MATRICES      ---------------------####
def initialise_matrix(shape):
    """Initialise numpy matrix of given shape with 0.

    Args:
        shape (tuple): shape of the matrix (nb of lines, nb of columns)

    Returns:
        np.ndarray: initialised matrix
    """
    return np.zeros(shape)


def fill_distance_matrix(template_prot):
    """Create a distance matrix for all pairwise distances in the template.

    Args:
        template_prot (Protein): template protein

    Returns:
        np.ndarray: distance matrix
    """
    n = template_prot.get_length()
    distance_matrix = initialise_matrix((n, n))
    # calculate all distances between alpha carbons in the protein
    for i in range(n):
        for j in range(n):
            if i != j:
                dist = template_prot.calculate_inter_ca_distance(i, j)
                distance_matrix[i, j] = round(dist, 2)
    return distance_matrix


def get_dope_score(dope_mat, res1, res2, distance):
    """Get corresponding dope score from dope dataframe.

    Args:
        dope_mat (pd.DataFrame): the dope scores dataframe
        res1 (str): name of the first residue (1-letter)
        res2 (str): name of the second residue
        distance (float): distance between the two residues in angstroms

    Returns:
        float: corresponding dope score
    """
    check_positive_number(distance)
    # we can't have only one atom - force infinity
    if distance == 0:
        return 0
    # return high value if distance too short or too high
    if distance < 0.25 or distance > 14.75:
        return 10
    # if we have the exact distance in the dope dataframe, return value
    if distance in dope_mat.columns:
        line_dope_mat = dope_mat[(dope_mat["res1"] == min(res1, res2)) &
                                 (dope_mat["res2"] == max(res1, res2))]
        return line_dope_mat[distance].iloc[0]
    # if we don't, apply an affine function
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
    score_inter = score_below + (((distance - distance_below) /
                                    (distance_above - distance_below)) *
                                    (score_above - score_below))
    return round(score_inter, 2)


def fill_LL_matrix(shape, dist_matrix, dope_matrix,
                   test_sequence, gap_penalty, i, j):
    """Fill and return a singular low-level matrix for the fixed point [i, j].

    Args:
        shape (tuple): shape of the matrix
        dist_matrix (np.ndarray): distance matrix for the template
        dope_matrix (pd.DataFrame): dataframe of the dope scores
        test_sequence (str): residues in the test sequence
        gap_penalty (int): gap penalty
        i (int): line number of the fixed point for this matrix
        j (int): column number of the fixed point for this matrix

    Returns:
        np.ndarray: low-level matrix for the fixed point [i, j]
    """
    # initialise matrix
    L_mat = initialise_matrix(shape)
    nrow = shape[0]
    ncol = shape[1]

    # initialise first line and column
    L_mat[1:, 0] = [k * gap_penalty for k in range(1, nrow)]
    L_mat[0, 1:] = [l * gap_penalty for l in range(1, ncol)]

    # go through matrix
    for k in range(1, nrow):
        for l in range(1, ncol):
            # if out of reach, add infinite number to force best path
            if (k < i and l >= j) or (k >= i and l < j):
                L_mat[k, l] = np.inf
            # skip lines already filled
            else:
                # calculate dope score
                dist = dist_matrix[j-1, l-1]
                res1 = test_sequence[i-1]
                res2 = test_sequence[k-1]
                dope_score = get_dope_score(dope_matrix, res1, res2, dist)
                # calculate minimum score
                L_mat[k, l] = min(
                    L_mat[k-1, l-1] + dope_score,
                    L_mat[k, l-1] + gap_penalty,
                    L_mat[k-1, l] + gap_penalty)
    #return matrix
    return L_mat


def create_global_LL_matrix(shape, dist_matrix, dope_matrix,
                            test_sequence, gap_penalty):
    """Create the 'super' matrix filled with each low-level matrix.

    Args:
        shape (tuple): shape of the matrix
        dist_matrix (np.ndarray): distance matrix for the template
        dope_matrix (pd.DataFrame): dataframe of the dope scores
        test_sequence (str): residues sequence in the test sequence
        gap_penalty (int): gap penalty

    Returns:
       np.ndarray: matrix with the low-level matrices in each of its cells
    """
    # initialise global matrix
    global_L_mat = np.empty(shape, dtype=object)
    # for each cell except first line and column, create a low-level matrix
    for i in range(1, shape[0]):
        for j in range(1, shape[1]):
            L_mat = fill_LL_matrix(shape, dist_matrix, dope_matrix,
                                   test_sequence, gap_penalty, i, j)
            global_L_mat[i, j] = L_mat
    return global_L_mat


def get_score_HL_matrix(global_L_mat, i, j):
    """Get the score of each low-level matrix to fill the high-level.

    Args:
        global_L_mat (np.ndarray): 'super' low-level matrix
        i (int): line number of the fixed point for the low-level matrix
        j (int): column number of the fixed point for the low-level matrix

    Returns:
        float: score of the corresponding low-level matrix
    """
    # for the corresponding low level matrix
    nrow = global_L_mat.shape[0]
    ncol = global_L_mat.shape[1]
    # return the LL matrix score
    return global_L_mat[i, j][nrow-1, ncol-1]


def fill_HL_matrix(shape, global_L_mat, gap_penalty):
    """Create and fill the high-level matrix, and store the path.

    Args:
        shape (tuple): shape of the matrix
        global_L_mat (np.ndarray): 'super' low-level matrix
        gap_penalty (int): gap penalty

    Returns:
        H_mat (np.ndarray): the high-level matrix
        align_mat (np.ndarray) : the matrix with one of the optimum paths
    """
    # initialise matrix
    H_mat = initialise_matrix(shape)
    align_mat = initialise_matrix(shape)
    nrow = shape[0]
    ncol = shape[1]
    
    # initialise first line and column
    H_mat[1:, 0] = [k * gap_penalty for k in range(1, nrow)]
    H_mat[0, 1:] = [l * gap_penalty for l in range(1, ncol)]

    # run through high level matrix
    for i in range(1, shape[0]):
        for j in range(1, shape[1]):
            # get the corresponding low-level score
            score = get_score_HL_matrix(global_L_mat, i, j)
            # calculate minimum score
            H_mat[i, j] = min(H_mat[i-1, j-1] + score,
                              H_mat[i, j-1] + gap_penalty,
                              H_mat[i-1, j] + gap_penalty)
            # mark corresponding path in the other matrix
            if H_mat[i, j] == H_mat[i-1, j-1] + score:
                align_mat[i, j] = 1
            elif H_mat[i, j] == H_mat[i, j-1] + gap_penalty:
                align_mat[i, j] = 2
            elif H_mat[i, j] == H_mat[i-1, j] + gap_penalty:
                align_mat[i, j] = 3
            else:
                print(f"Error: high-level matrix not valid")
                break
    return H_mat, align_mat


def backtracking(align_mat, sequence, template_prot):
    """Backtrack high-level matrix to get aligment.

    Args:
        align_mat (np.ndarray) : matrix with one of the optimum paths
        sequence (str): residues in the test sequence
        template_prot (Protein): template protein

    Returns:
        str: one of the optimum alignments
    """
    seq_aligned_reverse = ""
    res_aligned_reverse = ""
    # start from end point and go back in matrix
    i = align_mat.shape[0] - 1
    j = align_mat.shape[1] - 1
    while i > 0 and j > 0:
        # if diagonal, add match
        if align_mat[i, j] == 1:
            seq_aligned_reverse += sequence[i-1]
            res_aligned_reverse += str(j)[::-1]
            i -= 1
            j -= 1
        # if up, add gap to query
        elif align_mat[i, j] == 2:
            seq_aligned_reverse += "-"
            res_aligned_reverse += str(j)[::-1]
            j -= 1
        # if left, add gap to template
        elif align_mat[i, j] == 3:
            seq_aligned_reverse += sequence[i-1]
            res_aligned_reverse += "-"
            i -= 1
        else:
            print(f"Error: backtracking matrix non valid")
            break
    # return alignment in right order
    return f"{seq_aligned_reverse[::-1]}\n{res_aligned_reverse[::-1]}"


##############################################################################
#####################                MAIN                #####################
##############################################################################

if __name__ == "__main__":
    ####------------       IMPORT ARGS       ------------####
    args = read_args()
    seq_file = args.seq_file
    template_file = args.template_file
    GAP_PENALTY = args.gap_penalty

    ####------------       READ DATA       -------------####
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
    mat_shape = (len(row_names_query) + 1, n_atoms_template + 1)

    # build 'super' low-level-matrix
    global_L_mat = create_global_LL_matrix(mat_shape, dist_matrix,
                                           dope_scores, test_seq, GAP_PENALTY)
    # build high-level matrix
    H_mat, path_mat = fill_HL_matrix(mat_shape, global_L_mat, GAP_PENALTY)

    ####-----------    CALCULATE ALIGNMENT    -----------####
    alignment = backtracking(path_mat, test_seq, template_prot)

    ####--------------    PRINT RESULTS    --------------####
    print(alignment)
    print(f"Total energy : {round(H_mat[-1, -1], 2)}")
