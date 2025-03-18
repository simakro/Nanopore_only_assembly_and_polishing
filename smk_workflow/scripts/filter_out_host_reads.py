from time import localtime, strftime
from collections import defaultdict
import argparse
import subprocess


class Mapping:
    """class for storage and handling of mapping information"""

    def __init__(
        self,
        qname,
        qlen,
        qstart,
        qend,
        samestrand,
        tname,
        tlen,
        tstart,
        tend,
        matches,
        total_bp,
        qual,
        kwattr,
    ):
        self.qname = qname
        self.qlen = int(qlen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.samestrand = self.eval_strand(samestrand)
        self.tname = tname
        self.tlen = int(tlen)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.matches = int(matches)
        self.total_bp = int(total_bp)
        self.qual = int(qual)
        self.kwattr = kwattr
        self.gen_kw_attr()

    @staticmethod
    def eval_strand(strand_info):
        """evaluate strand info"""
        return True if strand_info == "+" else False

    def gen_kw_attr(self):
        """generate class attributes from key-worded entries in mapping output"""
        kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
        for key in kwattr_dict:
            self.__dict__[key] = kwattr_dict[key]


def get_args():
    parser = argparse.ArgumentParser()

    required = parser.add_argument_group("Required")
    required.add_argument(
        "-q", "--query", dest="query",
        help="Path to query fastaq file",
        required=True
    )
    required.add_argument(
        "-hgr",
        "--host_genome_ref",
        dest="host_genome_ref",
        help="Path to host genome reference file",
        required=True,
    )
    required.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Path to output file",
    )
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "-ft",
        "--filetype",
        dest="filetype",
        default="fastq",
        help="Type of nucl. acid sequence file (fasta/fastq) Default='fastq'",
    )
    args = parser.parse_args()
    return args


def minimap_subprocess(
    args: object,
    secondary: str="no",
    check: bool=True,
    capture_output: bool=True,
) -> str:
    out_paf = ".".join(args.query.split(".")[:-1]) + ".paf"
    cmd = [
        "minimap2",
        "-x",
        "map-ont",
        "-c",  # calculate cigar in order to get exact base matches number
        f"--secondary={secondary}",
        f"{args.host_genome_ref}",
        f"{args.query}",
        "-o", f"{out_paf}",
    ]
    subprocess.run(cmd, check=check, capture_output=capture_output)
    return out_paf


def load_paf(paf: str) -> dict:
    """load all paf data from file"""
    with open(paf, "r", encoding="utf-8") as paf:
        map_dct = defaultdict(list)
        for line in paf:
            if len(line.strip()) > 0:
                lisp = line.strip().split("\t")
                mapping = Mapping(*lisp[:12], lisp[12:])
                map_dct[mapping.qname].append(mapping)
    return map_dct


def get_mapped_readnames(map_dct: dict) -> list:
    mapped_reads = set()
    for rname in map_dct:
        mapped_reads.add(rname)
    print(f"Filtering out {len(mapped_reads)} reads mapping to host reference.")
    return list(mapped_reads)


class ReadExtractor:

    def __init__(
            self,
            args: object,
            read_names: list,
            mode: str = "inverse",
            filetype: str="fastq"
        ):
        self.mode = mode  # "normal" or "inverse"
        self.in_readfile: str = args.query
        self.outfile: str = args.output
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


def main():
    datestr = strftime("%Y-%m-%d_%H-%M-%S", localtime())
    print(f"Starting process filter_out_host_reads at {datestr}")
    print("Mapping reads against host genome to filter out host reads.")
    args: object = get_args()
    out_paf: str = minimap_subprocess(args)
    map_dct: dict = load_paf(out_paf)
    read_names: list = get_mapped_readnames(map_dct)
    ReadExtractor(args, read_names, mode="inverse", filetype=args.filetype)


if __name__ == "__main__":
    main()
