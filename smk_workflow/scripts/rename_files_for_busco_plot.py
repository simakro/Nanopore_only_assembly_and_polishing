import os
import sys
from glob import glob

def rename_busco_files(path):
    busco_sums = glob(os.path.join(path, "short_summary.specific.*"))
    for file in busco_sums:
        print(file)
        asmnm = "_".join(file.split("specific.")[1].split(".txt")[0].split("."))
        new_name =  f"short_summary.specific.bacteria.{asmnm}.txt"
        os.replace(file, os.path.join(path, new_name))

if __name__ == "__main__":
    path = sys.argv[1]
    rename_busco_files(path)
