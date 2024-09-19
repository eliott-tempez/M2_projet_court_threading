"""
This script is used to evaluate the output of meta_alignment.py.
It takes as input all text files in the directory passed as argument,
and returns for each of the protein of interest : a histogram of 
the total energy of the alignment to each template, and a histogram
of the gap proportion compared to the aligned residues. You can find the 
documentation at https://github.com/eliott-tempez/M2_projet_court_threading

Usage:
======
    python evaluate_alignment.py argument1
    
    argument1: the directory with outputs from the meta_alignment.py script

    Returns :
        a histogram of the total energy of the alignment to each template, 
        and a histogram of the gap proportion compared to the aligned residues
        for each protein of interest, automatically saved in the 
        results/eval/ folder.
        
"""

__author__ = "Eliott TEMPEZ"
__date__ = "2024-09-12"


import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


OUTPUT_DIR = "results/eval"


##############################################################################
###################                FUNCTIONS                ##################
##############################################################################

####-----------------------      HANDLE ERRORS      ----------------------####
def check_file(file_path):
    is_file = os.path.isfile(file_path)
    assert is_file, "Error : directory results/alignments should " \
    "only contain files"


####------------------------      READ ARGS      -------------------------####
def read_args():
    """Read and return command line arguments.

    Returns:
        argparse.Namespace: command lign arguments
    """
    descr = "Evaluate the output of meta_alignment.py"
    parser = argparse.ArgumentParser(description=descr)

    # arguments
    folder_align = "parent folder containing one or more txt files," \
        " generated from the program meta_alignment.py "

    parser.add_argument("folder_align", type=str, help=folder_align)
    args = parser.parse_args()
    return args


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


def calculate_gaps_proportion(code, alignment):
    """Calculate the proportion of gaps inserted in the alignment.

    Args:
        code (str): name of the protein
        alignment (str): alignment sequences

    Returns:
        float: gap proportion
    """
    # length of the sequence in the text file
    first_line_lst = alignment.split("\n")[0].split()
    len_seq = len(first_line_lst)
    second_line_lst = alignment.split("\n")[1].split()
    n_gaps = len_seq - len([i for i in second_line_lst if i == "|"])
    # calculate proportion
    return (n_gaps/len_seq)


def read_txt_file(file_path):
    """Read text file in parent folder and return a dictionnary.

    Args:
        file_path (str): text file path

    Returns:
        dict: dictionnary with all informations found in the file
    """
    protein_dict = {}
    with open(file_path, "r") as file_in:
        content = file_in.read()

        # extract each entry
        entries = content.strip().split("\n\n")
        for entry in entries:
            lines = entry.split("\n")
            # extract prot name and type
            line_0 = lines[0].split()
            if len(line_0) == 2:
                code, entry_type = lines[0].split()
            else:
                code, entry_type = lines[0].strip(), ""
            # initialise dict
            if code not in protein_dict:
                protein_dict[code] = {
                    "type": entry_type,
                    "alignment": "",
                    "score": 0.0,
                    "gap_prop": 0.0
                }

            # extract alignment
            alignment = "\n".join(lines[1:4])
            protein_dict[code]["alignment"] = alignment
            # extract proportion
            protein_dict[code]["gap_prop"] = calculate_gaps_proportion(code, alignment)
            # extract score
            score_line = lines[-1]
            score = float(score_line.split()[-1])
            protein_dict[code]["score"] = score
    return protein_dict


####--------------------------      PLOTS      ---------------------------####
def plot_histogram(protein_dict, prot_name, output_dir):
    """Plot histogram from each protein of interest and save the png file.

    Args:
        protein_dict (dict): dictionnary for the protein of interest
        prot_name (_type_): name of the protein of interest
        output_dir (_type_): output directory where to save the histogram file
    """
    templates = []
    scores = []
    props = []
    colors = []

    # define a color mapping for different types of proteins
    type_color_mapping = {
        "all_alpha": "skyblue",
        "all_beta": "salmon",
        "alpha+beta": "orchid",
        "other": "lightgrey"
    }
    # for each template protein
    for template_prot in protein_dict:
        if template_prot == prot_name:
            prot_type = protein_dict[template_prot]["type"]
        templates.append(template_prot)
        scores.append(protein_dict[template_prot]["score"])
        props.append(protein_dict[template_prot]["gap_prop"])
        # get the type of the protein and assign a color
        type = protein_dict[template_prot].get("type", "other")
        color = type_color_mapping.get(type, "lightgrey")
        colors.append(color)
    
    # for the global score
    # sort by scores (ascending order)
    sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i])
    templates = [templates[i] for i in sorted_indices]
    scores = [scores[i] for i in sorted_indices]
    colors = [colors[i] for i in sorted_indices]
    # create histogram
    plt.figure(figsize=(10, 6))
    plt.bar(templates, scores, color=colors)
    plt.xlabel("template protein name")
    plt.ylabel("score")
    plt.title(f"Alignment scores for sequence of interest {prot_name} ({prot_type})")
    plt.xticks(rotation=45)
    plt.tight_layout()
    # legend
    legend_handles = [mpatches.Patch(color=color, label=ptype) for ptype, color in type_color_mapping.items()]
    plt.legend(handles=legend_handles, title="Protein Type")
    # save the figure
    output_path = os.path.join(OUTPUT_DIR, f"histogram_{prot_name}.png")
    plt.savefig(output_path)
    plt.close()
    
    # for the gap insertions
    # sort by proportion (ascending order)
    sorted_indices = sorted(range(len(props)), key=lambda i: props[i])
    templates = [templates[i] for i in sorted_indices]
    props = [props[i] for i in sorted_indices]
    colors = [colors[i] for i in sorted_indices]
    # create histogram
    plt.figure(figsize=(10, 6))
    plt.bar(templates, props, color=colors)
    plt.xlabel("template protein name")
    plt.ylabel("gap proportion")
    plt.title(f"Gap proportion for sequence of interest {prot_name} ({prot_type})")
    plt.xticks(rotation=45)
    plt.tight_layout()
    # legend
    legend_handles = [mpatches.Patch(color=color, label=ptype) for ptype, color in type_color_mapping.items()]
    plt.legend(handles=legend_handles, title="Protein Type")
    # save the figure
    output_path = os.path.join(OUTPUT_DIR, f"histogram_prop_{prot_name}.png")
    plt.savefig(output_path)
    plt.close()


##############################################################################
#####################                MAIN                #####################
##############################################################################

if __name__ == "__main__":
    # import args
    args = read_args()
    INPUT_DIR = args.folder_align
    
    # get files list
    file_paths = extract_files_from_dir(INPUT_DIR)

    # get infos from files
    for file in file_paths:
        prot_name = get_prot_name_from_file(file)
        prot_dict = read_txt_file(file)
        # draw plot for each protein of interest
        plot_histogram(prot_dict, prot_name, OUTPUT_DIR)
