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


import pandas as pd
import numpy as np



DOPE_FILE = "data/dope.par"
#### A SUPPRIMER POUR AJOUTER LES COMMANDES UTILISATEURS
seq_file = "data/1UX8.fasta"
template_file = "data/1UX8.pdb"




def read_dope_file(filename):
    """Read dope file and return pandas dataframe.

    Args:
        filename (_str_): _dope file name_

    Returns:
        _pd.DataFrame_: _dataframe of dope score for each pairwise interaction_
    """
    # get the names of the distances corresponding to the dope scores
    anstrom_len_as_str = np.arange(0.75, 15.5, 0.5).tolist()
    # read dope file
    dope_mat_full = pd.read_csv(filename, sep=" ",
                                names=["res1", "atom1", "res2", "atom2"] + anstrom_len_as_str)
    # keep only carbon alpha dope scores
    dope_mat_ca = dope_mat_full[(dope_mat_full["atom1"] == "CA") & (dope_mat_full["atom2"] == "CA")]
    # return dataframe without unimportant columns
    return dope_mat_ca.drop(columns=["atom1", "atom2"]).reset_index()




if __name__ == "__main__":
    # read dope values
    dope_scores = read_dope_file(DOPE_FILE)
    
            

