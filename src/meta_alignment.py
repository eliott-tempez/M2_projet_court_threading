"""
This script is used to run the script align_structure.py on several 
templates, for one protein of interest. You can find the documentation at 
https://github.com/eliott-tempez/M2_projet_court_threading

Usage:
======
    python meta_alignment.py argument1 argument2

    argument1: a folder containing subfolders (optional), each containing
    one or more pdb files of the template proteins. If there are subfolder,
    their names should be some form of classification for the proteins
    argument2: the pdb name of the protein of interest. The corresponding
    fasta file must be in one of the subfolders mentioned above

    Returns :
        A text file, containing for each template protein:
        - The optimum alignment with the protein of interest
        - The total energy of this global alignment
"""

__author__ = "Eliott TEMPEZ"
__date__ = "2024-09-12"


OUTPUT_DIR = "results/"
# A RETIRER POUR INPUT LIGNE DE COMMANDE
folder = "data/poteins"
query = "1BAL"


import os



##############################################################################
###################                FUNCTIONS                ##################
##############################################################################

####-----------------------      HANDLE ERRORS      ----------------------####
def check_parentfold(has_subfolder, has_file, parent_dir):
    if has_subfolder and has_file:
        raise ValueError(f"{parent_dir} contains folder(s) and file(s) : \
            should contain only one of the two")
    elif not has_subfolder and not has_file:
        raise ValueError(f"{parent_dir} empty, should contain \
            folder(s) or file(s)")
        














####----------------------      HANDLE FOLDERS      ----------------------####
def check_subfolders(parent_dir):
    """Check if the directory from command line has subfolders.

    Args:
        parent_dir (str): general directory path

    Returns:
        bool: are there subfolders?
    """
    has_subfolder = False
    for element in os.listdir(parent_dir):
        # check if there's a subdirectory
        full_path = os.path.join(parent_dir, element)
        if os.path.isdir(full_path):
            has_subfolder = True
        else:
            has_file = True
    check_parentfold(has_subfolder, has_file, parent_dir)
    return has_subfolder