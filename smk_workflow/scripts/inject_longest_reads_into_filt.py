import os
import argparse
from get_x_cov_longest_reads import readfile_stats


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", "--acceptor",
        dest="acceptor",
        help="Path to acceptor file (typically the larger one)",
    )
    parser.add_argument(
        "-i", "--inject",
        dest="inject",
        help="Readsfile from which reads are going to be injected if they are not present in the acceptor file"
    )
    args = parser.parse_args()
    return args


def extract_reads(read_file, rn_lst, filetype, mode="normal"):
    header_ind = "@" if filetype == "fastq" else ">"
    outfile = read_file + ".reinject"
    with open(read_file, "r") as reads, open(outfile, "w") as extr:
        for line in reads:
            if line.startswith(header_ind):
                ls = line.strip().split(header_ind)[1].split(" ")
                rname = ls[0]
                if rname in rn_lst:
                    if mode == "inverse":
                        pass
                    else:
                        extr.write(line)
                        extr.write(next(reads))
                else:
                    if mode == "inverse":
                        extr.write(line)
                        extr.write(next(reads))
                    else:
                        pass
    return outfile

def reinject_long_reads(acceptor: str, select_reads: str):
    with open(acceptor, "a") as acc, open(select_reads, "r") as donor:
        inj = donor.read()
        acc.write(inj)
    root, ext = ".".join(acceptor.split(".")[:-1]), acceptor.split(".")[-1]
    mod_name = ".".join([root, "pluslong", ext])
    os.rename(acceptor, mod_name)
    return mod_name

def main():
    args = get_args()
    acceptor, inject = args.acceptor, args.inject
    # get dicts of readnames and lengths from fasta/qs
    acc_dict, acc_fformat, acc_log = readfile_stats(acceptor)
    inj_dict, inj_fformat, inj_log = readfile_stats(inject)
    # check which of the longest reads are already in the size and qual filtered reads
    intersect = set(acc_dict.keys()).intersection(set(inj_dict.keys()))
    print("Reads in acceptor file", len(acc_dict))
    print("Reads in injection file", len(inj_dict))
    print("Intersection of both files", len(intersect))
    reads_to_inject = [rn for rn in inj_dict.keys() if rn not in intersect]
    print("Reads to inject", len(reads_to_inject))
    # extract the long reads which have been filtered out by sq-filtering for reinjection
    selected_reads = extract_reads(inject, reads_to_inject, inj_fformat)
    mod_name = reinject_long_reads(acceptor, selected_reads)
    print(f"Renamed {acceptor} to {mod_name}")


if __name__ == "__main__":
    main()
