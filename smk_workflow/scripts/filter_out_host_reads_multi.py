from time import localtime, strftime, perf_counter
from collections import defaultdict
import multiprocessing as mp
import argparse
import subprocess
import os


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
        "-d", "--fq_dir", dest="fq_dir",
        help="Path to directory containing fastaq files",
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
        "--outdir",
        dest="outdir",
        help="Path to output directory",
    )
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "-ft",
        "--filetype",
        dest="filetype",
        default="fastq",
        help="Type of nucl. acid sequence file (fasta/fastq) Default='fastq'",
    )
    optional.add_argument(
        "-c",
        "--cores",
        dest="cores",
        type=int,
        default=1,
        help="Number of cores to be used as workers. Default=1",
    )
    
    args = parser.parse_args()
    return args


def minimap_subprocess(
    query: str,
    args: object,
    secondary: str="no",
    check: bool=True,
    capture_output: bool=True,
) -> str:
    qname = os.path.split(query)[1]
    paf_name = ".".join(qname.split(".")[:-1]) + ".paf"
    out_paf = os.path.join(args.outdir, paf_name)
    if not os.path.exists(out_paf):
        cmd = [
            "minimap2",
            "-x",
            "map-ont",
            "-c",  # calculate cigar in order to get exact base matches number
            f"--secondary={secondary}",
            f"{args.host_genome_ref}",
            f"{query}",
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


# def batch_extraction_proc(read_dir: str, args: object):
#     """
#     extract non-human reads from all individual fastq files in 
#     """
#     ilog.vlprint("Starting BLAST-searches for evaluation of int-sites", 1)

#     db_jobs = []
#     blj_info = {}
#     bl_jobs = []
#     error_log = open("sec_site_blastdb_err.log", "a")
#     for site in secondary_intsite_objs:
#         sec_site = secondary_intsite_objs[site]
#         db_path = os.path.join(args.outdir, "secsitedb_" + site + ".fa")
#         with open(db_path, "w") as db:
#             db.write(
#                 f">{sec_site.name} {sec_site.chrom}~{sec_site.spanned_region[0]}"
#                 f"-{sec_site.spanned_region[1]}\n"
#             )
#             db.write(sec_site.spanned_seq)
#         cif = InputFileChecker(db_path, "extensive", interact=False, silent=True)
#         corrupt = cif.corrupt
#         if not corrupt:
#             blj_info[site] = db_path
#             kw_args = {"msg": "secondary int_site candidate"}
#             dbp = mp.Process(target=gen_blastdb, args=(db_path, args), kwargs=kw_args)
#             db_jobs.append(dbp)
#             dbp.start()
#         else:
#             error_log.write(site + "\n")
#     error_log.close()

#     for proc in db_jobs:
#         proc.join()

#     for site in blj_info:
#         supp_reads = site + "_supporting_reads.fasta"
#         blp = mp.Process(
#             target=blast_eval_secondary, args=[site, blj_info[site], supp_reads, args]
#         )
#         bl_jobs.append(blp)
#         blp.start()

#     for proc in bl_jobs:
#         proc.join()

#     # get blast statistics with first site in blj_info
#     if len(list(blj_info.keys())) > 0:
#         test_sec_site = list(blj_info.keys())[0]
#         db_path = blj_info[test_sec_site]
#         stats_out = os.path.join(args.outdir, f"stats_{test_sec_site}_out")
#         blast_params = get_blast_stats(supp_reads, stats_out, db_path, args)
#     else:
#         blast_params = "not available"
#     # delete all blast-dbs in outdir-folder
#     cleanup_blastdb(args.outdir)
#     return blast_params



class ReadExtractor:

    def __init__(
            self,
            query: str,
            args: object,
            read_names: list,
            mode: str = "inverse",
            filetype: str="fastq"
        ):
        self.mode = mode  # "normal" or "inverse"
        self.in_readfile: str = query
        self.outdir: str = args.outdir
        self.outfile: str = self.get_outfile_name()  # args.outdir
        self.read_names: list = read_names
        self.filetype: str = filetype
        self.fastq = True if self.filetype == "fastq" else False
        self.header_ind = "@" if self.fastq else ">"
        self.l_ct = 0
        self.passed_ct = 0
        self.write_out_ct = 0
        self.reads = None
        self.extr = None
        if not os.path.exists(self.outfile):
            self.extract_nonhost_reads()

    def get_outfile_name(self):
        path = self.outdir
        in_name = os.path.split(self.in_readfile)[1]
        fname = ".".join(in_name.split(".")[:-1]) + "_nonhost.fastq"      
        fpath = os.path.join(path, fname)
        return fpath

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


def core(fastq, args):
    datestr = strftime("%Y-%m-%d_%H-%M-%S", localtime())
    print(f"Starting process filter_out_host_reads at {datestr}")
    print("Mapping reads against host genome to filter out host reads.")
    # args: object = get_args()
    out_paf: str = minimap_subprocess(fastq, args)
    map_dct: dict = load_paf(out_paf)
    read_names: list = get_mapped_readnames(map_dct)
    ReadExtractor(fastq, args, read_names, mode="inverse", filetype=args.filetype)


def concatenate_out_fastqs(args):
    # outdir = "results/{experiment}/{barcode}/nonhost_fastq",
    # outfile = "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq",
    cat_path: str = os.path.split(args.outdir)[0]
    res, exp, bc = cat_path.split(os.path.sep)
    cat_file: str = f"{exp}_{bc}_all_nonhost.fastq"
    cat_outfile: str = os.path.join(cat_path, cat_file)
    # cat_outfile = constr_outfile_path(args)
    fq2cat = [f.path for f in os.scandir(args.outdir) if f.path.endswith(".fastq")]
    if not os.path.exists(cat_outfile):
        with open(cat_outfile, "w") as out:
            for fq in fq2cat:
                with open(fq, "r") as fqo:
                    data = fqo.read()
                    out.write(data)
  

def main():
    args: object = get_args()
    # better use more flexible solution also allowing .fq etc.
    fastqs = [f.path for f in os.scandir(args.fq_dir) if f.path.endswith(".fastq")]
    os.makedirs(args.outdir, exist_ok=True)
    print(fastqs)
    # for fq in fastqs:
    #     core(fq, args) 
    with mp.Pool(processes=args.cores) as pool:
        arglist = [(fq, args) for fq in fastqs]
        pool.starmap(core, arglist)
    concatenate_out_fastqs(args)
    

if __name__ == "__main__":
    main()


    # class Args:
    #     def __init__(self, query_file, out_file):
    #         self.query = query_file
    #         self.outdir = out_file

    # args = Args(
    #     # "/home/simon/Nanopore_only_assembly_and_polishing/smk_workflow/data/Neuropatho_SM035_N473_plasmacytoma/N473_fastq_pass_all/FAR29577_pass_e3891f15_9a5312a6_5.fastq",
    #     "/home/simon/Nanopore_only_assembly_and_polishing/smk_workflow/results/Neuropatho_N1244_MB_Verdacht/fastq_pass/Neuropatho_N1244_MB_Verdacht_fastq_pass_1_all.fastq",
    #     "test_extraction.fastq"
    # )

    # read_names = [
    #     "a4db619c-c49f-42d3-a12e-d175c349f337",
    #     "09f62bdf-477f-4277-ace4-7ad5f722e6be",
    #     "2141da61-bc94-4c1e-a635-714fedd7b568",
    #     "b2750860-a5ec-4bc7-8193-a9eda8866dd5",
    #     "04eca36f-a985-40fd-8436-9351b7a352ad",
    #     "1c861362-6292-42b3-a527-c950cc00012e",
    #     "dc475a8d-f1e3-4fc2-8026-eb01a9f7f0ce",
    #     "a116eb99-ed8c-4c10-8e3d-9e72a31bd5b5",
    #     "7d1b1aa7-e4d2-468f-97c4-993d565d5d02",
    #     "9f53d70c-712f-4592-b0c7-7883adbe8be6",
    #     "25527a1d-785a-463f-8c65-47f282855aa1",
    #     "a41a786c-9af1-498b-9eb7-9a24f12b8b1e",
    #     "3f210677-c6b4-4857-a30e-fe13cb38d35b",
    #     "476e8169-7935-4fc5-8c5b-b82b93311565",
    #     "57517cb7-df11-453d-bf49-446c1c87389a",
    #     "8f822529-02e5-4c14-bf2a-58388ab2a5ad",
    #     "2472b7d9-38f9-4022-934d-3ae60a38b619",
    #     "768d122e-b3b1-4cad-9a88-02e75d632826",
    #     "42896e0f-5b2e-41b8-91eb-0884aa9c0bb2",
    #     "c992dbb0-b399-47c7-a865-fc794e48d72a",
    #     "bda368d4-0276-418c-8972-460596a4dee8",

    #     "95768083-54c2-4590-b957-4492205e5cad",
    #     "d34b9345-4fcf-4084-b3cc-265d006aa151",
    #     "f7270abe-329c-47d6-a6a8-1970f6b4cfb8",
    #     "618f4d7b-db53-4051-ba47-9a14769b692b"


    # ]

    # start = perf_counter()
    # for _ in range(1):
    #     # ReadExtractor(args, read_names, mode="inverse", filetype="fastq")
    #     ReadExtractor2(args, read_names, mode="inverse", filetype="fastq")
    # end = perf_counter()
    # runtime = end - start
    # print("runtime", runtime)



