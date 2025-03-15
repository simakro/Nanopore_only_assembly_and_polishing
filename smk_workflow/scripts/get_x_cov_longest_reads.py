import random
import argparse
import sys
import os
import statistics as stat
import math


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-cov",
        "--x_cov",
        dest="coverage",
        type=int,
        help="x-times coverage to be achieved with the longest available reads.",
    )
    parser.add_argument(
        "-g",
        "--genome_size",
        dest="genome_size",
        type=int,
        help="Size of genome in bases",
    )
    parser.add_argument(
        "-r", "--read_file", dest="read_file", help="Path to read-file (fasta/q)"
    )
    parser.add_argument(
        "-m",
        "--mode",
        default="normal",
        dest="mode",
        help="Extraction mode. Default=normal. In normal mode only the reads wi"
        "th the names in the list are selected. Type ‘inverse‘ to select all re"
        "ads, except the ones in the list.",
    )
    args = parser.parse_args()
    return args


def readfile_stats(readfile):
    # This function I took and from my own stats_and_split_reads_from_ec-centric.py with slight modifications
    # os.makedirs(out_dir, exist_ok=True)
    fastaq = readfile
    with open(fastaq, "r") as reads:
        read_dict = {}
        qstring_comp = {}
        readlen = []
        fformat = set()
        seq = ""
        header = ""
        lc = 0
        for line in reads:
            lc += 1
            # print(lc)
            if line.startswith("@"):
                fformat.add("fastq")
                header = line
                len_read = len(next(reads).strip())
                qstring_comp[lc] = set()
                qstring_comp[lc].add(len_read)
                readlen.append(len_read)
                read_dict[header] = len_read
                lc += 1
            elif line.startswith("+"):
                qstring_comp[lc - 2].add(len(next(reads).strip()))
                lc += 1
            elif line.startswith(">"):
                fformat.add("fasta")
                if lc > 1:
                    read_dict[header] = len(seq)
                    readlen.append(len(seq))
                seq = ""
                header = line
            else:
                if "fasta" in fformat:
                    seq = seq + line.strip()
        longest = max(readlen)
        print("longest", longest)
        longest_read = {k: v for (k, v) in read_dict.items() if v == longest}
        print("longest_read", longest_read)
        qstring_comp = {k: v for (k, v) in qstring_comp.items() if len(v) > 1}
        if len(fformat) == 1:
            if "fasta" in fformat:
                if lc > 2 * len(readlen) + 2:
                    print(
                        "WARNING: multiline fasta with est. line length:",
                        math.ceil(((sum(readlen) / (lc - len(readlen))) / 10)) * 10,
                    )
                    print(
                        "WARNING: get_x_longest_reads.py script is currently only optimized for single-line read file formats"
                    )
                    sys.exit()
                else:
                    print("single line fasta")
            else:
                print("Read file format = fastq")
        elif len(fformat) == 0:
            print("Could not detect any read file format")
            sys.exit()
        else:
            print(
                "WARNING: Multiple read file format characteristics detected:", fformat
            )
            print("WARNING: The chosen input file is likely corrupted.")
            print("WARNING: SM_readQC can not properly evaluate mixed files.")
            sys.exit()

        num_reads = len(readlen)
        # out_dir = os.path.split(readfile)[0]
        log_file = readfile + ".stats.txt"
        file_format = list(fformat)[0]
        with open(log_file, "w") as stat_log:
            print(f"File format: {file_format}")
            print(f"Number of lines in file: {lc}", file=stat_log)
            print(f"No. of reads {num_reads}", file=stat_log)
            print(f"Average read length: {stat.mean(readlen)}", file=stat_log)
            print(f"Median read length: {stat.median(readlen)}", file=stat_log)
            if "fasta" in fformat:
                print(f"longest read: {longest_read}", file=stat_log)
            if "fastq" in fformat:
                print(f"longest read: {max(readlen)}", file=stat_log)
                if len(qstring_comp) > 0:
                    print(
                        f"WARNING: The length of the quality string does not match"
                        f"readlength in {len(qstring_comp)} cases. {qstring_comp}",
                        file=stat_log,
                    )
            print(f"Total base count {sum(readlen)}", file=stat_log)

    return read_dict, file_format, log_file


def report_reads_to_extract(x_longest, len_ct, log_file=None):
    print(
        f"Selected {len(x_longest)} reads with a total of {len_ct} bases for extraction",
        file=log_file,
    )
    print("selected the following reads for extraction:", file=log_file)
    for r in x_longest:
        print(r[0].strip(), file=log_file)
        print(r[1], file=log_file)


def get_x_longest(read_dict, args, logfile):
    reads_by_len = sorted(list(read_dict.items()), key=lambda x: x[1])[::-1]
    required_bases = args.coverage * args.genome_size
    x_longest = []
    len_ct = 0
    for r in reads_by_len:
        if len_ct < required_bases:
            x_longest.append(r)
            len_ct += r[1]
        else:
            break
    with open(logfile, "a") as logfile:
        report_reads_to_extract(x_longest, len_ct, log_file=None)
        report_reads_to_extract(x_longest, len_ct, log_file=logfile)
    rn_lst = [rt[0] for rt in x_longest]
    return rn_lst


def parse_read_file(read_file, rn_lst, filetype, mode="normal"):
    header_ind = "@" if filetype == "fastq" else ">"
    with open(read_file, "r") as reads, open(read_file + ".longestx", "w") as extr:
        for line in reads:
            if line.startswith(header_ind):
                if line in rn_lst:
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


def main():
    args = get_args()
    read_dict, file_format, log_file = readfile_stats(args.read_file)
    read_name_lst = get_x_longest(read_dict, args, log_file)
    parse_read_file(args.read_file, read_name_lst, file_format, mode=args.mode)


if __name__ == "__main__":
    main()
