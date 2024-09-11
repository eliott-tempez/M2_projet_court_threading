"""
This script is used to run the script align_structure.py on several 
templates, for one protein of interest. You can find the documentation at 
https://github.com/eliott-tempez/M2_projet_court_threading

Usage:
======
    python meta_alignment.py argument1 argument2

    argument1: a parent folder containing one or more pdb files of 
    the template proteins. If there are subfolders (optional),
    their names should be some form of classification for the proteins
    argument2: path to the fasta file of the protein of interest (can be 
    in the parent folder or not)

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
def check_parentfold(has_subfolder, has_file, parent_path):
    """Check if general directory contains folder/file (not both or none)"""
    if has_subfolder and has_file:
        raise ValueError(f"{parent_path} contains folder(s) and file(s) : \
            should contain only one of the two")
    elif not has_subfolder and not has_file:
        raise ValueError(f"{parent_path} empty, should contain \
            folder(s) or file(s)")


def check_fasta_or_pdb(filepath):
    """Check if a file is a fasta/pdb file."""
    is_fasta_or_pdb = filepath.split(";")[-1] in ["fasta", "pdb"]
    assert is_fasta_or_pdb, "Error : file in (sub)folder should only \
        be fasta or pdb"

        














####----------------------      HANDLE FOLDERS      ----------------------####
def check_if_subfolders(parent_path):
    """Check directory path has subfolders.

    Args:
        parent_path (str): general directory path

    Returns:
        bool: are there subfolders?
    """
    has_subfolder = False
    # check if parent path is a folder
    if not os.path.isdir(parent_path):
        raise ValueError("First argument should be path to a directory")
    # check each element in parent path
    for element in os.listdir(parent_path):
        # check if there's a subdirectory
        full_path = os.path.join(parent_path, element)
        if os.path.isdir(full_path):
            has_subfolder = True
        else:
            has_file = True
    check_parentfold(has_subfolder, has_file, parent_path)
    return has_subfolder


def list_subfolders(parent_path):
    """Read and list subfolders in parent directory.

    Args:
        parent_path (str): general directory path containing subfolders

    Returns:
        list: list of the names of the subfolders
    """
    subfolders = []
    for subfold in os.listdir(parent_path):
        subfolders.append(subfold)
        sub_path = os.path.join(parent_path, subfold)
        for element in sub_path:
            full_path = os.path.join(sub_path, element)
            check_fasta_or_pdb(full_path)
            if os.path.isdir(full_path):
                raise ValueError("Subfolders should not contain folders")
    return subfolders


def extract_pdb_files_from_dir(dir_path):
    """Extract pdb filenames from a directory.

    Args:
        dir_path (str): directory path

    Returns:
        list: list of the pdb filenames
    """
    pdb_filenames = []
    for file in os.listdir(dir_path):
        if file.split(".")[-1] == "pdb":
            pdb_filenames.append(file)
    return pdb_filenames


def get_all_pdb_paths(parent_path):
    """Get all paths to pdb files in parent folder.

    Args:
        parent_path (str): general directory path

    Returns:
        list: list of all pdb paths
    """
    pdb_paths = []
    has_subfolders = check_if_subfolders(parent_path)
    if has_subfolders:
        subfolders = list_subfolders(parent_path)
        for subfold in subfolders:
            sub_path = os.path.join(parent_path, subfold)
            pdb_paths += [os.path.join(sub_path, pdbfile) for 
                          pdbfile in extract_pdb_files_from_dir(sub_path)]
    else:
        pdb_paths += [os.path.join(parent_path, pdbfile) for
                      pdbfile in extract_pdb_files_from_dir(parent_path)]
    return pdb_paths
        















if __name__ == "__main__":
    # get all pdb paths
    pdb_paths = get_all_pdb_paths(folder)
    # if subfolders, get their names
    if check_if_subfolders(folder):
        subs = list_subfolders(folder)