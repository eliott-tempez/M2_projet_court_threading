"""
This script is used to evaluate the output of meta_alignment.py.
You can find the documentation at 
https://github.com/eliott-tempez/M2_projet_court_threading

Usage:
======
    python evaluate_alignment.py

    Returns :
        A set of graphs, automatically saved in the results/eval directory
"""

__author__ = "Eliott TEMPEZ"
__date__ = "2024-09-12"


import os
import matplotlib.pyplot as plt


INPUT_DIR = "results/alignments/"
OUTPUT_DIR = "results/eval"


##############################################################################
###################                FUNCTIONS                ##################
##############################################################################

####-----------------------      HANDLE ERRORS      ----------------------####
def check_file(file_path):
    is_file = os.path.isfile(file_path)
    assert is_file, "Error : directory results/alignments should " \
    "only contain files"




####-----------------------      HANDLE FILES      -----------------------####
def extract_files_from_dir(dir_path):
    """Extract file paths from a directory.

    Args:
        dir_path (str): directory path

    Returns:
        list: list of the file paths
    """
    file_paths = []
    for file in os.listdir(dir_path):
        file_path = os.path.join(dir_path, file)
        check_file(file_path)
        file_paths.append(file_path)
    return file_paths


def get_prot_name_from_file(file_path):
    """Get the protein name from a file path.

    Args:
        file_path (str): file path

    Returns:
        str: protein name (PDB format)
    """
    file_name = file_path.split("/")[-1]
    sans_txt = file_name.split(".")[0]
    prot_name = sans_txt.split("_")[-1]
    return prot_name












##############################################################################
#####################                MAIN                #####################
##############################################################################

if __name__ == "__main__":
    # get files list
    file_paths = extract_files_from_dir(INPUT_DIR)