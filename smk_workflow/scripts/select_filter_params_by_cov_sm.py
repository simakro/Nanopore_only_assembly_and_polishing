# from Nanopore_only_assembly_and_polishing
# this script/tool would benefit from a better fastq-file integrity check
# as is, fastq files lacking + and qstring are not recognized as corrupted,
# but judged to be empty

import os
import sys
from time import perf_counter as timer
# from timeit import timeit
import statistics as stat
import json
import argparse
from math import log


def get_args():
    parser = argparse.ArgumentParser()

    required = parser.add_argument_group("Required")
    required.add_argument("-f", "--fastq", help="path to fastq read file", type=str, required=True)
    required.add_argument("-s", "--genome_size", help="Size of target genome", type=int, required=True)

    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "-c", "--cov",
        default=30, help="Specify coverage aimed for [default=30]", type=int
    )
    optional.add_argument(
        "-q", "--qual",
        default=10, help="Request read quality aimed for [default=10]", type=float
    )
    optional.add_argument(
        "-l", "--len",
        default=10000, help="Request read length aimed for [default=10000]", type=int
    )
    # optional.add_argument(
    #     "-o", "--out_dir",
    #     default="", help="Define directory for results json (& dumped reads if set) file [default=current working directory]", type=str
    # )
    optional.add_argument(
        "-d", "--dump",
        default=False, choices=["fasta", "fastq"],
        help="Set this flag to dump the selected reads in a file of chosen format. [default=False; choices=(fastq/fasta)]",
        # action="store_true"
    )
    # optional.add_argument(
    #     "-ofmt", "--out_format",
    #     default="fastq", help="Set output format for dumped reads. [default=fastq]", type=str
    # )
    optional.add_argument(
        "-opt", "--optimize", type=typechk_bool("optimize"), default=True,
        help="Toggle optimization of read quality and length cutoffs. If optimization is deactivated, coverage will be ignored. [default=True; choices=(True/False)]", 
    )
    optional.add_argument(
        "-decall", "--allow_reduction_all",
        default=False, help="Allow to reduce requested read length and quality parameters until requested coverage is reached [default=False]", 
        action="store_true"
    )
    # optional.add_argument(
    #     "-declen", "--allow_reduction_len",
    #     default=False, help="Allow to reduce requested read length parameter until requested coverage is reached [default=False]", 
    #     action="store_true"
    # )
    # optional.add_argument(
    #     "-decqual", "--allow_reduction_qual",
    #     default=False, help="Allow to reduce requested read quality parameter until requested coverage is reached [default=False]", 
    #     action="store_true"
    # )
    optional.add_argument(
        "-lss", "--len_stepsize",
        default=1000, help="During read-selection optimization use this stepsize for length [default=1000]", 
    )
    optional.add_argument(
        "-qss", "--qual_stepsize",
        default=1, help="During read-selection optimization use this stepsize for quality [default=1]", 
    )

    silent = parser.add_argument_group("Silent")
    silent.add_argument(
        "--direction",
        default=None, help=argparse.SUPPRESS 
    )

    args = parser.parse_args()
    return args


# def parse_and_analyze_fastq_org(fastq):
#     fstart = timer()
#     bases = ["A", "T", "G", "C"]
#     ascii_chars =  "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
#     qchars = set(ascii_chars)
#     qchar_dict = {str(q): ord(q)-33 for q in qchars}

#     used_qchars = set()
#     av_read_quality = []
#     start_baseq = {num+1:list() for num in range(100)}
#     # print(start_baseq)
#     end_baseq = {-(num+1):list() for num in range(100)}
#     # print(end_baseq)

#     # fastq must be in single-line format for this to work
#     with open(fastq, "r") as fq:
#         header = False
#         seq_line = False
#         qstring = False
#         for line in fq:
#             if line.startswith("@") and not qstring:
#                 header = True
#                 qstring = False
#             elif line == ("+\n") and seq_line:
#                 seq_line = False
#                 qstring = True
#             elif line[0].upper() in bases and header: # fastq must be in single-line format for this to work
#                 header = False
#                 seq_line = True
#             elif line[0] in qchars and qstring:
#                 qstring = False # fastq must be in single-line format for this to work
#                 lset = set(line.strip())
#                 used_qchars = used_qchars.union(lset)
#                 lq = [ord(c)-33 for c in line.strip()]
#                 for n in range(len(lq[:100])):
#                     start_baseq[n+1].append(lq[n])
#                 for n in range(len(lq[:100])):
#                     end_baseq[-(n+1)].append(lq[-(n+1)])
#                 mean_qs = stat.mean(lq)
#                 av_read_quality.append(mean_qs)

#     av_quality_reads = stat.mean(av_read_quality)
#     sd_qual_reads = stat.stdev(av_read_quality)


#     print(f"Average read quality = {av_quality_reads}")
#     print(f"SD read quality = {sd_qual_reads}")
#     for n in start_baseq:
#         start_baseq[n] = stat.mean(start_baseq[n])
#     print(f"Average start_baseq for all reads: {start_baseq}")
#     for n in end_baseq:
#         end_baseq[n] = stat.mean(end_baseq[n])
#     print(f"Average end_baseq for all reads: {end_baseq}")
#     fend = timer()
#     print(f"Runtime {parse_and_analyze_fastq_org.__name__}: {fend-fstart}")
#     # print(sorted(qchar_dict.items()))

#     return  start_baseq, end_baseq


class Read():
    def __init__(self, header: str, len: int, av_qual: float):
        self.header = header
        self.qual = av_qual
        self.len = len


def typechk_bool(arg_name):
    """Type checker for boolean cli arguments"""
    
    def chk_bool_arg(arg_val):
        try:
            arg = eval(arg_val)
        except TypeError as te:
            print(te)
            raise argparse.ArgumentTypeError(
                f"Accepted values for {arg_name} are 'False' and 'True'."
                )
        except NameError as ne:
            print(ne)
            raise argparse.ArgumentTypeError(
                f"Accepted values for {arg_name} are 'False' and 'True'."
                )
        
        if type(arg)!=bool:
            print(f"{arg_val} is no valid value for argument {arg_name}")
            raise argparse.ArgumentTypeError(
                f"Accepted values for {arg_name} are 'False' and 'True'."
                )
        else:
            return arg
    
    return chk_bool_arg


def parse_fastq(fastq: str, analyze=True, search_lst=False, outfile=False, outfmt=False) -> list:
    fstart = timer()
    # bases = ["A", "T", "G", "C", "N"]
    bases = ["A", "T", "G", "C", "N", "a", "t", "g", "c", "n"]
    ascii_chars =  "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    qchars = set(ascii_chars)
    qchar_dict = {str(qc): ord(qc)-33 for qc in qchars}
    err_prob = {qc:10**(-qs/10) for qc,qs in qchar_dict.items()}

    reads = []

    if search_lst:
        if not outfile:
            print(
                "No outfile defined. If parse_fastq is used in search mode, outfile must be provided."
                )
            sys.exit()
        else:
            search_lst = [r.header for r in search_lst]
            out = open(outfile, "w")

    with open(fastq, "r") as fq:
        curr_header = None
        curr_seqline = None
        header = False
        seq_line = False
        qstring = False
        dump = False
        for line in fq:
            if line.startswith("@") and not qstring:
                curr_header = line
                header = True
                qstring = False
                if search_lst:
                    if curr_header in search_lst:
                        dump = True
                    else:
                        dump = False
            elif line == ("+\n") and seq_line:
                seq_line = False
                qstring = True
            elif line[0] in bases and header:
                header = False
                seq_line = True
                if dump:
                    curr_seqline = line
            elif line[0] in qchars and qstring:
                line = line.strip()
                qstring = False
                if analyze:
                    lq = [err_prob[q] for q in line.strip()]
                    mean_qs = -10 * log(sum(lq)/len(lq), 10)
                    reads.append(Read(curr_header, len(line), mean_qs))
                if dump:
                    if outfmt=="fastq":
                        out.write(curr_header)
                        out.write(curr_seqline)
                        out.write("+\n")
                        out.write(line + "\n")
                    else:
                        out.write(">" + curr_header[1:])
                        out.write(curr_seqline)
                    dump = False

    if search_lst:
        out.close()               

    fend = timer()
    print(f"Runtime {parse_fastq.__name__}: {fend-fstart}")
    if analyze:
        if len(reads) > 0:
            print("Input read file statistics:")
            read_ct = len(reads)
            read_lens = [read.len for read in reads]
            total_bp = sum(read_lens)
            print(f"Read count: {read_ct}")
            print(f"Total bases: {total_bp}")
            # if read_ct > 0:
            print(f"Avg read length: {total_bp/read_ct}")
            # else:
            #     print(f"Avg read length: 0")
            print(f"Shortest read: {min(read_lens)}")
            print(f"Longest read: {max(read_lens)}")
        else:
            print("No reads were available")
        # ReadN50, median length ...
    return  reads


def optimize_cutoffs(    
    org_reads: list,
    selection: list,
    target_cov: float,
    genome_size: int,
    req_qual: float,
    req_len: int,
    args: object,
    direction=None):
    selection = [read for read in org_reads if read.len>=req_len and read.qual>=req_qual]
    total_bp = sum([read.len for read in selection])
    read_num = len(selection)
    total_cov = round(total_bp/genome_size, 1)
    print(f"total_cov: {total_cov}, read_num: {read_num}, qual: {req_qual}, len: {req_len}")
    if total_cov < target_cov:
        if len(selection) < len(org_reads):
            if not args.direction=="up":
                if args.allow_reduction_all:
                    args.direction = "down"
                    print("More reads available for selection, decrementing len and qual parameters.")
                    print(req_qual-args.qual_stepsize, req_len-args.len_stepsize)
                    req_qual, req_len, selection = optimize_cutoffs(org_reads, selection, target_cov, genome_size, req_qual-args.qual_stepsize, req_len-args.len_stepsize, args)
                    return req_qual, req_len, selection
                else:
                    print("Target coverage could not be reached with the set qual and len limits. If you want to allow downward adjustment of these parameters use the -decall/--allow_reduction_all flag.")
                    return req_qual, req_len, selection
            else:
                print("Target coverage threshold was crossed during parameter optimization. Returning last good parameter set.")
                req_qual, req_len = req_qual-args.qual_stepsize, req_len-args.len_stepsize
                selection = [read for read in org_reads if read.len>=req_len and read.qual>=req_qual]
                return req_qual, req_len, selection
        else:
            print("All reads have been selected, but target coverage could not be reached")
            print(req_qual, req_len)
            return req_qual, req_len, selection
    elif total_cov > target_cov:
        if args.direction=="down":
            print(f"Target coverage could be reached by decreasing selection parameters to qual: {req_qual}, len: {req_len}")
            return req_qual, req_len, selection
        else:
            args.direction = "up"
            # incr_qual = optimize_selection(reads, target_cov, genome_size, req_qual+args.qual_increment, req_size)
            # incr_len = optimize_selection(reads, target_cov, genome_size, req_qual, req_size+args.len_increment)
            print(f"More coverage than requested ({total_cov}/{target_cov}), incrementing len and qual parameters.")
            req_qual, req_len, selection = optimize_cutoffs(org_reads, selection, target_cov, genome_size, req_qual+args.qual_stepsize, req_len+args.len_stepsize, args)
            # trials = [incr_qual, incr_len, incr_both]
            return req_qual, req_len, selection
    else:
        print("Exactly reached target coverage")
        return req_qual, req_len, selection


def select_reads(
    reads: list,
    target_cov: float,
    genome_size: int,
    req_qual: float,
    req_len: int,
    args: object,
    direction=None):
    print("Running selection optimization")
    fstart = timer()
    # if req_len<=1000 and args.len_stepsize > 100:
    #     args.len_stepsize = 100
    selection = [read for read in reads if read.len>=req_len and read.qual>=req_qual]
    total_bp = sum([read.len for read in selection])
    read_num = len(selection)
    total_cov = round(total_bp/genome_size, 1)
    print(f"total_cov: {total_cov}, read_num: {read_num}, qual: {req_qual}, len: {req_len}")
    if args.optimize:
        req_qual, req_len, selection = optimize_cutoffs(reads, selection, target_cov, genome_size, req_qual, req_len, args)
        total_bp = sum([read.len for read in selection])
        read_num = len(selection)
        total_cov = round(total_bp/genome_size, 1)
        print("Parameter optimization for  cutoffs completed.")
        print(f"Returning reads with qual >= {req_qual} and length >= {req_len}.")
        print(f"total_bp: {total_bp} (~{total_cov}x cov), read_count: {read_num}.")
    else:
        print(f"Returning all reads with qual >= {req_qual} and length >= {req_len} without optimization of cutoffs.")
        print(f"total_bp: {total_bp} (~{total_cov}x cov), read_count: {read_num}.")
    fend = timer()
    print(f"Runtime {select_reads.__name__}: {fend-fstart}")
    return req_qual, req_len, selection


def run_cover_up():
    fstart = timer()
    args = get_args()
    reads = parse_fastq(args.fastq)
    rqual, rlen, selection = select_reads(
        reads,
        args.cov,
        args.genome_size,
        args.qual,
        args.len,
        args
        )
    print(f"Final values: rqual={rqual}, rlen={rlen}")
    out_file = os.path.join(os.path.split(args.fastq)[0], "filt_params.json")
    with open(out_file, "w") as out:
        data = json.dumps({"len": rlen, "qual": rqual})
        out.write(data)
    if args.dump:
        file_name = ".".join(os.path.split(args.fastq)[1].split(".")[:-1]) + f"_sqfilt." + args.dump
        reads_out = os.path.join(os.path.split(args.fastq)[0], file_name)
        parse_fastq(args.fastq, analyze=False, search_lst=selection, outfile=reads_out, outfmt=args.dump)
    fend = timer()
    print(f"Total runtime: {fend-fstart}")
    return rqual, rlen

if __name__ == "__main__":
    rqual, rlen = run_cover_up()
