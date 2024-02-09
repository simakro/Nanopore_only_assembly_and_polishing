import os
import sys
import json
from glob import glob
import shutil

"""
Template for requested dict/json very general:
##############################################
{
    "/path/to/source_dir": {
        "folders": {
            "folder_name": {
                    "files": [],
                    "folders": {},
            },
        },
        "files": [],
    },
}

Template for requested dict/json more verbose:
#################################
{
    "/path/to/source_dir": {
        "folders": {
            "specific_folder_name": {
                "folders": {
                    "Data1": {
                        "files": ["*"],
                        "folders": {}
                    },
                },
                "files":[],
            },
            "recurring_pattern*": {
                "folders": {
                    "Data1": {
                        "files": ["spec_file1", spec_file2],
                        "folders": {"*"}
                    },
        "files": [],
    },
}

Theoretical example for requested dict/json:
################################
{
    "/home/user/repo/workflow/results/experiment": {
        "folders": {
            "barcode*": {
                "folders": {
                    "tool1": {
                        "files": ["*.png"],
                    },
                    "tool2": {
                        "folders": {
                            "tool2_data": {
                                "files": ["*.tsv", "*.json", "*.txt"],
                            },
                        },
                    },
                    "tool3": {
                        "folders": {
                            "tool3_data": {
                                "files": ["medaka_flye_pilon3.oriented.fasta"],
                            },
                        },
                    },
                    "tool4": {"files": ["*"]},
                    "tool5": {"files": ["*.html"]},
                    "tool6": {"files": ["*.svg", "*.png"]},
                },
            },
        },
        "files": [],
    },
}
"""

# requested = {
#     "/home/simon/Nanopore_only_assembly_and_polishing/smk_workflow/results/SM0026_bacterial_genomes_KP_curr": {
#         "folders": {
#             "barcode*": {
#                 "folders": {
#                     "assemblytics": {
#                         "files": ["*.png"],
#                     },
#                     "busco_medaka_flye_pilon3": {
#                         "folders": {
#                             "run_bacteria_odb10": {
#                                 "files": ["*.tsv", "*.json", "*.txt"],
#                             },
#                         },
#                     },
#                     "circl_fixstart": {
#                         "folders": {
#                             "medaka_flye_pilon3": {
#                                 "files": ["medaka_flye_pilon3.oriented.fasta"],
#                             },
#                         },
#                     },
#                     "dnadiff": {"files": ["*"]},
#                     "fastqc_nanopore": {"files": ["*.html"]},
#                     "ilmn_clip_sm": {"files": ["*.html"]},
#                     "ilmn_fastc_before": {"files": ["*.html"]},
#                     "medaka_flye_pilon3_circos": {"files": ["*.svg", "*.png"]},
#                     "medaka_flye_pilon3_gtdbtk": {"files": ["*.tsv", "*.json", "*.txt"]},
#                     "medaka_flye_pilon3_plasclass": {"files": ["*.out"]},
#                     "medaka_flye_pilon3_update-sinfo": {"files": ["*.json"]},
#                     "prokka_medaka_flye_protgbk_class": {"files": ["*.gff", "*.tbl", "*.tsv", "*.txt"]},
#                     "nanoplot": {"files": ["*report.html"]},
#                     "flye":{"files": ["*info.txt", "*.gv"]}, # , ".gfa"
#                 },
#             },
#             "busco_graph": {
#                 "folders": {
#                     "*": {"files": ["*.png"]}
#                 },
#             },
#             "medaka_flye_pilon3_gtdbtk_sinfo": {
#                 "files": ["*json"]
#             },
#         },
#         "files": ["filt_params.json"],
#     },
# }


requested = {
    "/home/simon/Nanopore_only_assembly_and_polishing/smk_workflow/results/SM0026_bacterial_genomes_KP": {
        "folders": {
            "barcode*": {
                "folders": {
                    "flye":{"files": ["*info.txt", "*.gv", "*.gfa"]}, # , "*.gfa"
                },
            },
        },
    },
}


def walk_req_tree(requested):
    # print("requested", requested)
    final_paths = []
    for source_folder in requested:
        growing_path = [source_folder]
        construct_paths(
            requested[source_folder],
            growing_path,
            final_paths
            )
    # print ("Ultimate paths", final_paths)
    return final_paths


def construct_paths(source_folder, growing_path, final_paths):
    for entry in source_folder:
        if entry=="files":
            curr_path = os.path.join(*growing_path)
            for file in source_folder["files"]:
                final_paths.append(os.path.join(curr_path, file))
        elif entry=="folders":
            curr_path = os.path.join(*growing_path)
            for folder in source_folder["folders"]:
                curr_branch = list(growing_path)
                curr_branch.append(folder)
                final_paths = construct_paths(source_folder["folders"][folder], curr_branch, final_paths)
    return final_paths


def load_requested_json(req_js_path: str):
    with open(req_js_path, "r") as jsn_in:
        data = json.loads(jsn_in.read())
    return data


def expand_paths(constr_paths):
    expanded_paths = []
    for path in constr_paths:
        glob_paths = glob(path)
        expanded_paths.extend(glob_paths)
    print("expanded_paths", expanded_paths)
    return expanded_paths


def copy_to_summary_dir(expanded_paths, requested, summary_dir):
    for source_dir in requested:
        expath = [f for f in expanded_paths if source_dir in f]
        for source_path in expath:
            sub_path = source_path.split(source_dir)[1]
            target_path = "".join([summary_dir, sub_path])
            os.makedirs(os.path.split(target_path)[0], exist_ok=True)
            shutil.copy(source_path, target_path)


if __name__ == "__main__":
    summary_dir = sys.argv[1]
    constr_paths = walk_req_tree(requested)
    exp_paths = expand_paths(constr_paths)
    copy_to_summary_dir(exp_paths, requested, summary_dir)
