import os
import sys
from time import perf_counter as timer
from timeit import timeit
import statistics as stat
import json
import argparse


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
    optional.add_argument(
        "-o", "--out_dir",
        default="", help="Define directory for results json file [default=current working directory]", type=str
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



def parse_and_analyze_fastq_org(fastq):
    fstart = timer()
    bases = ["A", "T", "G", "C"]
    ascii_chars =  "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    qchars = set(ascii_chars)
    qchar_dict = {str(q): ord(q)-33 for q in qchars}

    used_qchars = set()
    av_read_quality = []
    start_baseq = {num+1:list() for num in range(100)}
    # print(start_baseq)
    end_baseq = {-(num+1):list() for num in range(100)}
    # print(end_baseq)

    # fastq must be in single-line format for this to work
    with open(fastq, "r") as fq:
        header = False
        seq_line = False
        qstring = False
        for line in fq:
            if line.startswith("@") and not qstring:
                header = True
                qstring = False
            elif line == ("+\n") and seq_line:
                seq_line = False
                qstring = True
            elif line[0].upper() in bases and header: # fastq must be in single-line format for this to work
                header = False
                seq_line = True
            elif line[0] in qchars and qstring:
                qstring = False # fastq must be in single-line format for this to work
                lset = set(line.strip())
                used_qchars = used_qchars.union(lset)
                lq = [ord(c)-33 for c in line.strip()]
                for n in range(len(lq[:100])):
                    start_baseq[n+1].append(lq[n])
                for n in range(len(lq[:100])):
                    end_baseq[-(n+1)].append(lq[-(n+1)])
                mean_qs = stat.mean(lq)
                av_read_quality.append(mean_qs)
                
    av_quality_reads = stat.mean(av_read_quality)
    sd_qual_reads = stat.stdev(av_read_quality)


    print(f"Average read quality = {av_quality_reads}")
    print(f"SD read quality = {sd_qual_reads}")
    for n in start_baseq:
        start_baseq[n] = stat.mean(start_baseq[n])
    print(f"Average start_baseq for all reads: {start_baseq}")
    for n in end_baseq:
        end_baseq[n] = stat.mean(end_baseq[n])
    print(f"Average end_baseq for all reads: {end_baseq}")
    fend = timer()
    print(f"Runtime {parse_and_analyze_fastq_org.__name__}: {fend-fstart}")
    # print(sorted(qchar_dict.items()))

    return  start_baseq, end_baseq


class Read():
    def __init__(self, len: int, av_qual: float):
        self.qual = av_qual
        self.len = len


def parse_and_analyze_fastq(fastq):
    fstart = timer()
    bases = ["A", "T", "G", "C"]
    ascii_chars =  "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    qchars = set(ascii_chars)
    qchar_dict = {str(q): ord(q)-33 for q in qchars}

    reads = []

    # fastq must be in single-line format for this to work
    with open(fastq, "r") as fq:
        header = False
        seq_line = False
        qstring = False
        for line in fq:
            if line.startswith("@") and not qstring:
                header = True
                qstring = False
            elif line == ("+\n") and seq_line:
                seq_line = False
                qstring = True
            elif line[0].upper() in bases and header: # fastq must be in single-line format for this to work
                header = False
                seq_line = True
            elif line[0] in qchars and qstring:
                qstring = False # fastq must be in single-line format for this to work
                lq = [qchar_dict[c] for c in line.strip()]
                mean_qs = sum(lq)/len(lq)
                reads.append(Read(len(lq), mean_qs))

    fend = timer()
    print(f"Runtime {parse_and_analyze_fastq.__name__}: {fend-fstart}")

    return  reads

def optimize_selection(reads: list, target_cov: float, genome_size: int, req_qual: float, req_len: int, args: object, direction=None):
    print("Running selection optimization")
    fstart = timer()
    # if req_len<=1000 and args.len_stepsize > 100:
    #     args.len_stepsize = 100
    selection = [read for read in reads if read.len>=req_len and read.qual>=req_qual]
    total_bp = sum([read.len for read in selection])
    total_cov = total_bp/genome_size
    print(f"total_cov {total_cov}, qual: {req_qual}, len: {req_len}")
    fend = timer()
    print(f"Runtime {optimize_selection.__name__}: {fend-fstart}")
    if total_cov < target_cov:
        if len(selection) < len(reads):
            if not args.direction=="up":
                if args.allow_reduction_all:
                    args.direction = "down"
                    print("More reads available for selection, decrementing len and qual parameters.")
                    print(req_qual-args.qual_stepsize, req_len-args.len_stepsize)
                    req_qual, req_len = optimize_selection(reads, target_cov, genome_size, req_qual-args.qual_stepsize, req_len-args.len_stepsize, args)
                    return req_qual, req_len
                else:
                    print("Target coverage could not be reached with the set qual and len limits. If you want to allow downward adjustment of these parameters use the -decall/--allow_reduction_all flag.")
                    return req_qual, req_len
            else:
                print("Target coverage threshold was crossed during parameter optimization. Returning last good parameter set.")
                return req_qual-args.qual_stepsize, req_len-args.len_stepsize
        else:
            print("All reads have been selected, but target coverage could not be reached")
            print(req_qual, req_len)
            return req_qual, req_len
    elif total_cov > target_cov:
        if args.direction=="down":
            print(f"Target coverage could be reached by decreasing selection parameters to qual: {req_qual}, len: {req_len}")
            return req_qual, req_len
        else:
            args.direction = "up"
            # incr_qual = optimize_selection(reads, target_cov, genome_size, req_qual+args.qual_increment, req_size)
            # incr_len = optimize_selection(reads, target_cov, genome_size, req_qual, req_size+args.len_increment)
            print(f"More coverage than requested {total_cov}/{target_cov}, incrementing len and qual parameters.")
            req_qual, req_len = optimize_selection(reads, target_cov, genome_size, req_qual+args.qual_stepsize, req_len+args.len_stepsize, args)
            # trials = [incr_qual, incr_len, incr_both]
            return req_qual, req_len
    else:
        print("Exactly reached target coverage")
        return req_qual, req_len


def run_cover_up():
    args = get_args()
    reads = parse_and_analyze_fastq(args.fastq)
    rqual, rlen = optimize_selection(reads, args.cov, args.genome_size, args.qual, args.len, args)
    print(f"Final values: rqual={rqual}, rlen={rlen}")
    out_file = os.path.join(os.path.split(args.fastq)[0], "filt_params.json")
    with open(out_file, "w") as out:
        data = json.dumps({"len": rlen, "qual": rqual})
        out.write(data)
    return rqual, rlen

if __name__ == "__main__":
    rqual, rlen = run_cover_up()
    
