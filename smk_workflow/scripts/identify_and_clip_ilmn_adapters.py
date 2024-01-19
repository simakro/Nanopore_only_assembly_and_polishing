import os
import sys
import gzip
import argparse
from time import perf_counter
from collections import Counter, defaultdict
import multiprocessing as mp


def get_args():
    parser = argparse.ArgumentParser()
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument(
        "-r", "--read_files", required=True, nargs="+",
        help="Path to a single or two "
        "paired fastq file/s (--r {file1} {file2}). Gzipped accepted. No "
        "whitespaces in file paths/names allowed."
        )
    required_args.add_argument(
        "-a", "--adapters", required=True, help="Path to adapters fasta file"
        )
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument(
        "-o", "--outdir", default=False, help="Specify target directory for "
        " output files. Will be created, if not yet existing."
        )
    optional_args.add_argument(
        "-m", "--min_len", default=20, type=int, help="Mininmum length for "
        "clipped reads. Shorter reads will be excluded. [20]"
        )
    optional_args.add_argument(
        "-cs", "--clip_start", default=0, type=int,
        help="Cut n bases from the beginning of each read. [0]"
        )
    optional_args.add_argument(
        "-sh", "--show_histo", default=False, action="store_true",
        help="Set this flag to output histograms of adapter positions in reads."
        )
    optional_args.add_argument(
        "-pd", "--probing_depth", default=250000, type=int,
        help="Define how many reads from the beginning of each read fastq are "
        "probed with all available adapters in order to determine the adapter "
        "set that is used to scan the whole file. [250000]"
        )
    args = parser.parse_args()
    return args


def extract_gz(gzipped_file, keep_gz=True):
    unzipped_file = ".".join(gzipped_file.split(".")[:-1])
    with open(unzipped_file, "w") as extracted:
                    gio = gzip.open(gzipped_file)
                    for line in gio:
                        extracted.write(line.decode("utf-8"))
    if not keep_gz:
        os.remove(gzipped_file)
    return unzipped_file


def load_adapter(adapter_fa):
    adapter_dict = dict()
    with open(adapter_fa, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                header = line.strip()
                seq = next(fa).strip().upper()
                adapter_dict[seq] = header
    return adapter_dict


def sl_fastq_generator(fastq):
    with open(fastq, "r") as fq:
        for line in fq:
            if line.startswith("@"):
                yield (
                    line.strip(),
                    next(fq).strip().upper(),
                    next(fq).strip(),
                    next(fq).strip(),
                    )


# def evaluate_read_files(read_arg):
#     read_files = read_arg.split(" ")
#     print(read_files)
def run_parallel_clipping(args):
    return_dct = mp.Manager().dict()
    read_fq = args.read_files
    clip_jobs = []
    ct = 0
    for fq in read_fq:
        ct += 1
        clj = mp.Process(
            target=chk_clip_ilmn_reads,
            args=(fq, args, return_dct),
            # kwargs=kw_args
            )
        clip_jobs.append(clj)
        print(f"Starting job {ct}")
        clj.start()
    for proc in clip_jobs:
        proc.join()
    print("return_dct", return_dct)
    return return_dct


def probe_for_adapters(read_fq, args):
    # further speed-up would be achievable if in cases where multiple adapters are found,
    # if one or more of those adapters are always preceeded by a different adapter,
    # it would suffice to only return the preceeding adapter, since this would lead
    # to removal of the following ones as well. A shorter list for searching of each
    # seq line leads to significant speed-up. 
    adapters_count = Counter()
    adapter_dict = load_adapter(args.adapters)
    fa_gen = sl_fastq_generator(read_fq)
    read_ct = 0
    for rid,seq,plus,qual in fa_gen:
        read_ct += 1
        if read_ct > args.probing_depth:
            del(fa_gen)
            break
        for ad in adapter_dict:
            if ad in seq:
                adapters_count[ad] += 1
    print(f"Probed {os.path.split(read_fq)[1]} for adapters:\n{adapters_count}")
    return {k:adapter_dict[k] for k in adapters_count}


def chk_clip_ilmn_reads(read_fq, args, return_dct):
    # set arguments
    # read_fq = args.read_files
    ## adapter_fa = args.adapters
    clip_start= args.clip_start
    min_len = args.min_len
    show_ad_histo = args.show_histo
    # initialize containers, counters and switches
    ## adapter_dict = load_adapter(adapter_fa)
    read_dict = dict()
    ad_count = Counter()
    ad_pos = defaultdict(Counter)
    in_read_len = Counter()
    out_read_len = Counter()
    extracted = False
    ad_clipped_reads = {}
    excluded_lst = []
    read_ct = 0
    # decompress if necessary
    if read_fq.endswith(".gz"):
        print("decompressing")
        read_fq = extract_gz(read_fq, keep_gz=True)
        extracted = read_fq
    # set up fastq generator
    print("Setting up generator")
    adapter_dict = probe_for_adapters(read_fq, args)
    fa_gen = sl_fastq_generator(read_fq)
    # define outfiles
    # clipped_out = ".".join(read_fq.split(".")[:-1]) + "_clipped.fq"
    # excluded_out = ".".join(read_fq.split(".")[:-1]) + "_excluded.fq"
    if args.outdir:
        print("Outdir arg is set")
        if not os.path.exists(args.outdir):
            print("Outdir does not exist")
            os.makedirs(args.outdir, mode=777)
            os.chmod(args.outdir, 0o777)
        else:
            print("Outdir already exists")
        path, name = os.path.split(read_fq)
        clip_name = ".".join(name.split(".")[:-1]) + "_clipped.fq"
        excl_name = ".".join(name.split(".")[:-1]) + "_excluded.fq"
        clipped_out = os.path.join(args.outdir, clip_name)
        excluded_out = os.path.join(args.outdir, excl_name)
    else:
        clipped_out = ".".join(read_fq.split(".")[:-1]) + "_clipped.fq"
        excluded_out = ".".join(read_fq.split(".")[:-1]) + "_excluded.fq"
    # if args.gzip_output:
    #     pass
    # iterate over fq entries in context manager
    with open(clipped_out, "w") as out, open(excluded_out, "w") as excl:
        print("Starting iteration over reads")
        for rid,seq,plus,qual in fa_gen:
            # set in-loop counters, switches and containers
            read_ct += 1
            adap_present = False
            ad_clipped = False
            exclude = False
            curr_adapters = {}
            # check if adapters are present in read
            for ad in adapter_dict:
                if ad in seq:
                    adap_present = True
                    ad_count[ad] += 1
                    pos = seq.index(ad)
                    ad_pos[ad][pos] += 1
                    curr_adapters[ad] = pos
            orig_len = len(seq)
            in_read_len[orig_len] += 1 # len of input reads; not reported yet
            # if clipping required, clip reads
            if any([adap_present, clip_start]):
                if adap_present:
                    min_pos = min(curr_adapters.items(), key=lambda x: x[1])[1]
                    ad_clipped = True
                else:
                    min_pos = orig_len
                seq, qual = seq[clip_start:min_pos], qual[clip_start:min_pos]
            # document clipping characteristics on read level
            if ad_clipped:
                clipped_len = len(seq)
                ad_clipped_reads[rid] = {
                    "adapter": curr_adapters,
                    "clip_start": clip_start,
                    "clipped_len": clipped_len,
                    "orig_len": orig_len,
                    "ad_bases_clipped": orig_len - min_pos,
                    "total_bases_clipped": orig_len - clipped_len,
                    }
                if clipped_len < min_len:
                    exclude = True
                    excluded_lst.append(rid)
            # write out read data to outfiles
            read_data = [rid,seq,plus,qual]
            if exclude:
                for item in read_data:
                    excl.write(f"{item}\n")
            else:
                for item in read_data:
                    out.write(f"{item}\n")
    # report stats and results
    for ad in ad_count:
        print(
            f"{adapter_dict[ad]} ({ad}) is present {ad_count[ad]} times in"
            f" read-fq"
            )
        if show_ad_histo:
            for n in range(max(ad_pos[ad])+1):
                col = int((ad_pos[ad][n]/20))
                if col > 0:
                    print(f"{n}\t{'*' * col}")
    
    ct_adclipped = len(ad_clipped_reads)
    bases_removed_adclip = [
        ad_clipped_reads[r]["ad_bases_clipped"] for r in ad_clipped_reads
        ]
    # all_removed_adclip = total_bases removed during adapter clipping, 
    # including everything beyond the adapter (e.g. aditional adapters/indices 
    # or sequencing artifacts like PolyG)
    all_removed_adclip = sum(bases_removed_adclip) 
    print(f"Processed a total of {read_ct} reads.")
    print(
        f"Cut adapters from {ct_adclipped} removing a total of "
        f"{all_removed_adclip} bases."
        )
    print(
        f"Excluded {len(excluded_lst)} reads with post-clip length < {min_len}"
        )
    if ct_adclipped > 0:
        avg_removed_adclip = all_removed_adclip/ct_adclipped
        print(
            f"On average removed {avg_removed_adclip} bases per read during"
            f" adapter clipping."
            )
    if clip_start:
        read_start_clip = clip_start*read_ct
        print(
            f"Cut {clip_start} from start of {read_ct} reads removing a total "
            f"of {read_start_clip} bases.")
        print(f"Removed a total of {all_removed_adclip + read_start_clip}")
    # remove unpacked fq if unzipping took place before (keeping only orig. gz)
    if extracted:
        os.remove(extracted)
    # return excluded_lst
    return_dct[clipped_out] = excluded_lst


def compare_excluded(return_dct, args):
    # if args.interleaved:
    #     pass
    # if args.not_paired:
    #     pass
    # else:
    fqs1 = set([r.split(" ")[0] for r in return_dct.values()[0]])
    fqs2 = set([r.split(" ")[0] for r in return_dct.values()[1]])
    asym = fqs1.symmetric_difference(fqs2)
    excl_asym = {
        k:[
            r.split(" ")[0] for r in v if r.split(" ")[0] in asym
            ] for k,v in return_dct.items()
        }
    # switch the exclusion lists
    if len(excl_asym)==2:
        key_lst = list(excl_asym.keys())
        switch_dct = dict(zip(key_lst, key_lst[::-1]))
        tmp_dict = {switch_dct[k]:v for k,v in excl_asym.items()}
        excl_asym = tmp_dict 
    else:
        print(
            "WARNING: excl_asym unexpectedly does not contain exactly 2 "\
            "file-keys as expected for paired-end non-interleaved fastqs. "\
            "Asymmetries in fastqs, if present, can not be fixed.‚"
            )

    return excl_asym


def repair_files(clipped_fq, excl_lst, args):
    print(f"Checking {clipped_fq}. Excl.-List: {excl_lst}")
    repaired = clipped_fq + ".rep"
    excl_fq = "_".join(clipped_fq.split("_")[:-1]) + "_excluded.fq"
    fa_gen = sl_fastq_generator(clipped_fq)
    with open(clipped_fq, "r") as clip, open(repaired, "w") as rep, open(excl_fq, "a") as excl:
        for rid,seq,plus,qual in fa_gen:
            read_data = [rid,seq,plus,qual]
            if rid.split(" ")[0] in excl_lst:
                print(
                    f"Removed {rid} from {clipped_fq} to restore pe-fq symmetry"
                    )
                for item in read_data:
                    excl.write(f"{item}\n")
            else:
                for item in read_data:
                    rep.write(f"{item}\n")
    os.remove(clipped_fq)
    os.replace(repaired, clipped_fq)


if __name__ == "__main__":
    args = get_args()
    start = perf_counter()
    return_dct = run_parallel_clipping(args)
    excl_asym = compare_excluded(return_dct, args)
    for fq in excl_asym:
        repair_files(fq, excl_asym[fq], args)
    # excluded_lst = chk_clip_ilmn_reads(args)
    stop = perf_counter()
    print(f"Runtime: {stop-start}")