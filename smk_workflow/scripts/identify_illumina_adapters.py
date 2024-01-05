import os
import sys
from time import perf_counter
from collections import Counter


def load_adapter(adapter_fa):
    adapter_dict = dict()
    with open(adapter_fa, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                header = line.strip()
                seq = next(fa).strip()
                adapter_dict[seq] = header
    return adapter_dict


def sl_fastq_generator(fasta):
    with open(fasta, "r") as fa:
        for line in fa:
            if line.startswith("@"):
                yield (line.strip(), next(fa).strip())


def chk_ilmn_reads(read_fa, adapter_fa):
    adapter_dict = load_adapter(adapter_fa)
    read_dict = dict()
    ad_count = Counter()
    fa_gen = sl_fastq_generator(read_fa)
    for id,seq in fa_gen:
        for ad in adapter_dict:
            if ad in seq:
                # read_dict[id] = seq
                ad_count[ad] += 1
    for ad in ad_count:
        print(f"{adapter_dict[ad]} ({ad}) is present {ad_count[ad]} times in read-fq")
    return ad_count


if __name__ == "__main__":
    adapter_fa = sys.argv[1]
    read_fa = sys.argv[2]
    start = perf_counter()
    ad_count = chk_ilmn_reads (read_fa, adapter_fa)
    stop = perf_counter()
    print(f"Runtime: {stop-start}")