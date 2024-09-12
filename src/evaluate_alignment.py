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
import matplotlib.patches as mpatches


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
                        "score": 0.0
                    }

            # extract alignment
            alignment = "\n".join(lines[1:4])
            protein_dict[code]["alignment"] = alignment
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
        templates.append(template_prot)
        scores.append(protein_dict[template_prot]["score"])
        # get the type of the protein and assign a color
        type = protein_dict[template_prot].get("type", "other")
        color = type_color_mapping.get(type, "lightgrey")
        colors.append(color)
        
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
    plt.title(f"Alignment scores for sequence of interest {prot_name}")
    plt.xticks(rotation=45)
    plt.tight_layout()
    # legend
    legend_handles = [mpatches.Patch(color=color, label=ptype) for ptype, color in type_color_mapping.items()]
    plt.legend(handles=legend_handles, title="Protein Type")
    # save the figure
    output_path = os.path.join(OUTPUT_DIR, f"histogram_{prot_name}.png")
    plt.savefig(output_path)
    plt.close()


##############################################################################
#####################                MAIN                #####################
##############################################################################

if __name__ == "__main__":
    # get files list
    file_paths = extract_files_from_dir(INPUT_DIR)

    # get infos from files
    for file in file_paths:
        prot_name = get_prot_name_from_file(file)
        prot_dict = read_txt_file(file)
        # draw plot for each protein of interest
        plot_histogram(prot_dict, prot_name, OUTPUT_DIR)
