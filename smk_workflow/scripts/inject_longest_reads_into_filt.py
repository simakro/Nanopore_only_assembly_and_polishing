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


class ReadExtractor:

    def __init__(
            self,
            args: object,
            read_names: list,
            mode: str = "inverse",
            filetype: str="fastq"  # "fasta"
            ):
        self.mode = mode  # "normal" or "inverse"
        self.in_readfile: str = args.inject
        self.outfile: str = args.inject + ".inject"
        self.read_names: list = read_names
        self.filetype: str = filetype
        self.fastq = True if self.filetype == "fastq" else False
        self.header_ind = "@" if self.fastq else ">"
        self.l_ct = 0
        self.passed_ct = 0
        self.write_out_ct = 0
        self.reads = None
        self.extr = None
        self.extract_nonhost_reads()

    def write_out_reads(self, line):
        self.write_out_ct += 1
        self.l_ct += 1
        self.extr.write(line)
        self.extr.write(next(self.reads))
        if self.fastq:
            self.l_ct += 2
            self.extr.write(next(self.reads))
            self.extr.write(next(self.reads))

    def extract_nonhost_reads(self):
        with open(self.in_readfile, "r") as self.reads, open(
            self.outfile, "w"
        ) as self.extr:
            for line in self.reads:
                self.l_ct += 1
                if line.startswith(self.header_ind) and self.l_ct % 2 == 1:
                    identifier = line.strip().split(" ")[0][1:]
                    if identifier in self.read_names:
                        if self.mode == "inverse":
                            self.passed_ct += 1
                        else:
                            self.write_out_reads(line)
                    else:
                        if self.mode == "inverse":
                            self.write_out_reads(line)
                        else:
                            self.passed_ct += 1
        written_out = self.write_out_ct if self.mode == "normal" else self.passed_ct
        print(
            f"Wrote {written_out} reads of the total of {self.write_out_ct + self.passed_ct} reads to outfile"
        )


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
    selected_reads = ReadExtractor(
        args, reads_to_inject, mode="normal", filetype="inj_fformat"
    ).outfile
    # selected_reads = extract_reads(inject, reads_to_inject, inj_fformat)
    mod_name = reinject_long_reads(acceptor, selected_reads)
    print(f"Renamed {acceptor} to {mod_name}")


if __name__ == "__main__":
    main()
