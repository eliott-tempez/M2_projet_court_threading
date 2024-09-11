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


import argparse
import os
import subprocess


OUTPUT_DIR = "results/"


##############################################################################
###################                FUNCTIONS                ##################
##############################################################################

####-----------------------      HANDLE ERRORS      ----------------------####
def check_parentfold(has_subfolder, has_file, parent_path):
    """Check if general directory contains folder/file (not both or none)"""
    if has_subfolder and has_file:
        raise ValueError(f"{parent_path} contains folder(s) and file(s) : " \
            "should contain only one of the two")
    elif not has_subfolder and not has_file:
        raise ValueError(f"{parent_path} empty, should contain " \
            "folder(s) or file(s)")


def check_fasta_or_pdb(filepath):
    """Check if a file is a fasta/pdb file."""
    is_fasta_or_pdb = filepath.split(".")[-1] in ["fasta", "fa", "pdb"]
    assert is_fasta_or_pdb, "Error : file in (sub)folder should only " \
        "be fasta or pdb"


####------------------------      READ ARGS      -------------------------####
def read_args():
    """Read and return command line arguments.

    Returns:
        argparse.Namespace: command lign arguments
    """
    descr = "Run the script align_structure.py on several " \
        "templates, for one protein of interest"
    parser = argparse.ArgumentParser(description=descr)

    # arguments
    folder_descr = "parent folder containing one or more pdb files. If there" \
        " are subfolders, their names should be some form of " \
        "classification for the model proteins"
    fasta_descr = "fasta file path of the protein sequence you want to model"

    parser.add_argument("main_folder", type=str, help=folder_descr)
    parser.add_argument("fasta_path", type=str, help=fasta_descr)
    args = parser.parse_args()
    return args


####-----------------------      HANDLE FILES      -----------------------####
def check_if_subfolders(parent_path):
    """Check directory path has subfolders.

    Args:
        parent_path (str): general directory path

    Returns:
        bool: are there subfolders?
    """
    has_subfolder = False
    has_file = False
    # check if parent path is a folder
    if not os.path.isdir(parent_path):
        raise ValueError("First argument should be path to a directory")
    # check each element in parent path
    for element in os.listdir(parent_path):
        # if not hidden file
        if not element.startswith("."):
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
        for element in os.listdir(sub_path):
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
    pdb_paths = {}
    has_subfolders = check_if_subfolders(parent_path)
    if has_subfolders:
        subfolders = list_subfolders(parent_path)
        for subfold in subfolders:
            sub_path = os.path.join(parent_path, subfold)
            pdb_paths[subfold] = [os.path.join(sub_path, pdbfile) for
                                  pdbfile in extract_pdb_files_from_dir(sub_path)]
    else:
        pdb_paths["parent"] = [os.path.join(parent_path, pdbfile) for
                               pdbfile in extract_pdb_files_from_dir(parent_path)]
    return pdb_paths


def get_protein_name(file_path):
    """Get the protein name from a file path.

    Args:
        file_path (str): file path

    Returns:
        str: protein name (PDB format)
    """
    file_name = file_path.split("/")[-1]
    prot_name = file_name.split(".")[0]
    return prot_name
   

####---------------------      RUN MAIN PROGRAM      ---------------------####
def align_structures(pdb_paths, fasta_path):
    query_name = get_protein_name(fasta_path)
    # run align_structure.py for all templates
    with open(OUTPUT_DIR + f"alignment_{query_name}.txt", "w") as f_out:
        f_out.write(f"{query_name}\n\n")
        for key in pdb_paths:
            for pdb_path in pdb_paths[key]:
                template_name = get_protein_name(pdb_path)
                command = ["python", "src/align_structure.py", fasta_path, pdb_path]
                result = subprocess.run(command, capture_output=True, text=True)
                if key != "parent":
                    f_out.write(f"{template_name} {key}\n")
                else:
                    f_out.write(f"{template_name}\n")
                f_out.write(f"{result.stdout}\n")
             


##############################################################################
#####################                MAIN                #####################
##############################################################################

if __name__ == "__main__":
    # import args
    args = read_args()
    main_folder = args.main_folder
    fasta_path = args.fasta_path
    # get all pdb paths
    pdb_paths = get_all_pdb_paths(main_folder)
    # run program
    align_structures(pdb_paths, fasta_path)
