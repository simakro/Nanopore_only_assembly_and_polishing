from time import localtime, strftime
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


def write_out_reads(line, l_ct, inf_handle, outf_handle, write_out_ct, fastq):
    write_out_ct += 1
    l_ct += 1
    outf_handle.write(line)
    outf_handle.write(next(inf_handle))
    if fastq:
        l_ct += 2
        outf_handle.write(next(inf_handle))
        outf_handle.write(next(inf_handle))


def extract_reads(read_file, rn_lst, filetype, mode="normal"):
    fastq = True if filetype == "fastq" else False
    header_ind = "@" if filetype == "fastq" else ">"
    outfile = read_file + ".reinject"
    l_ct = 0
    passed_ct = 0
    write_out_ct = 0
    with open(read_file, "r") as reads, open(outfile, "w") as extr:
        for line in reads:
            l_ct += 1
            if line.startswith(header_ind) and l_ct % 2 == 1:
                ls = line.strip().split(header_ind)[1].split(" ")
                rname = ls[0]
                if rname in rn_lst:
                    if mode == "inverse":
                        passed_ct += 1
                    else:
                        write_out_reads(line, l_ct, reads, extr, write_out_ct, fastq)
                        # extr.write(line)
                        # extr.write(next(reads))
                else:
                    if mode == "inverse":
                        write_out_reads(line, l_ct, reads, extr, write_out_ct, fastq)
                        # extr.write(line)
                        # extr.write(next(reads))
                    else:
                        passed_ct += 1
    print(
        f"Extracted {write_out_ct} reads of the total of previously selected"
        f" {write_out_ct + passed_ct} reads for reinjection."
    )
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
    datestr = strftime("%Y-%m-%d_%H-%M-%S", localtime())
    print(f"Starting process inject_longest_reads_into_filt at {datestr}")
    print("Reinject previously extracted longest reads into quality filtered reads.")
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
