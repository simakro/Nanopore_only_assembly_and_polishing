import os
import sys
import csv
import unicodedata
from glob import glob
from io import StringIO


def rename_busco_files(path):
    busco_sums = glob(os.path.join(path, "short_summary.specific.*"))
    for file in busco_sums:
        asmnm = "_".join(file.split("specific.")[1].split(".txt")[0].split("."))
        asmnm = analyze_asm_name(asmnm)
        new_name =  f"short_summary.specific.bacteria.{asmnm}.txt"
        os.replace(file, os.path.join(path, new_name))


def load_bact_genera():
    with open("resources/List_bacterial_genera.csv", "r", encoding="utf-8-sig") as in_csv:
        data = in_csv.read()
        data_norm = StringIO(unicodedata.normalize("NFKC", data))
        dct_csv = list(csv.DictReader(data_norm, delimiter=";"))
        gen_set = set([d['Genus'] for d in dct_csv])
        return gen_set


def analyze_asm_name(asm_nm):
    bc_synonyms = ["barcode", "bc"]
    bact_gen = load_bact_genera()
    if any([syn in asm_nm for syn in bc_synonyms]):
        return asm_nm
    elif any([segm in bact_gen for segm in asm_nm.split("_")]):
        splt_lst = asm_nm.split("_")
        genus = [segm for segm in splt_lst if segm in bact_gen]
        gen_idx = splt_lst.index(genus[0])
        new_name = f"{splt_lst[gen_idx]}_{splt_lst[gen_idx+1]}"
        return new_name
    else:
        return asm_nm


if __name__ == "__main__":
    path = sys.argv[1]
    rename_busco_files(path)
