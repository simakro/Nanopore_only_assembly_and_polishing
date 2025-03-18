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
    # args.logfile = get_log_name(args)
    return args


# def get_log_name(args: object) -> str:
#     pass


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
    # with open(read_lst, "w", encoding="utf-8") as rel:
    for rname in map_dct:
        mapped_reads.add(rname)
    # print("mapped_reads", mapped_reads)
    print(f"Filtering out {len(mapped_reads)} reads mapping to host reference.")
    return list(mapped_reads)


def write_out_reads(line, l_ct, inf_handle, outf_handle, write_out_ct, fastq):
    write_out_ct += 1
    l_ct += 1
    outf_handle.write(line)
    outf_handle.write(next(inf_handle))
    if fastq:
        l_ct += 2
        outf_handle.write(next(inf_handle))
        outf_handle.write(next(inf_handle))


def extract_nonhost_reads(args, rn_lst, mode="inverse"):
    fastq = True if args.filetype == "fastq" else False
    header_ind = "@" if fastq else ">"
    l_ct = 0
    passed_ct = 0
    write_out_ct = 0
    with open(args.query, "r") as reads, open(args.output, "w") as extr:
        for line in reads:
            l_ct += 1
            if line.startswith(header_ind) and l_ct % 2 == 1:
                identifier = line.strip().split(" ")[0][1:]
                if identifier in rn_lst:
                    if mode == "inverse":
                        passed_ct += 1
                    else:
                        write_out_reads(line, l_ct, reads, extr, write_out_ct, fastq)
                        # write_out_ct += 1
                        # extr.write(line)
                        # extr.write(next(reads))
                        # if fastq:
                        #     extr.write(next(reads))
                        #     extr.write(next(reads))
                else:
                    if mode == "inverse":
                        write_out_reads(line, l_ct, reads, extr, write_out_ct, fastq)
                        # write_out_ct += 1
                        # extr.write(line)
                        # extr.write(next(reads))
                        # if fastq:
                        #     extr.write(next(reads))
                        #     extr.write(next(reads))
                    else:
                        passed_ct += 1
    print(f"Wrote {write_out_ct} reads of the total of {write_out_ct + passed_ct} reads to outfile")


def main():
    datestr = strftime("%Y-%m-%d_%H-%M-%S", localtime())
    print(f"Starting process filter_out_host_reads at {datestr}")
    print("Mapping reads against host genome to filter out host reads.")
    args: object = get_args()
    out_paf: str = minimap_subprocess(args)
    map_dct: dict = load_paf(out_paf)
    read_names: list = get_mapped_readnames(map_dct)
    extract_nonhost_reads(args, read_names, mode="inverse")


if __name__ == "__main__":
    main()
