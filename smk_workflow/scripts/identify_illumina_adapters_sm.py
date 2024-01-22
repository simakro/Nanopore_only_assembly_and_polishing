import os
import sys
import gzip
from time import perf_counter
from collections import Counter, defaultdict


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
    # ad_pos = defaultdict(set)
    ad_pos = defaultdict(Counter)
    extracted = False
    if read_fa.endswith(".gz"):
        read_fa = extract_gz(gzipped_file, keep_gz=True)
        extracted = read_fa
    fa_gen = sl_fastq_generator(read_fa)
    for id,seq in fa_gen:
        for ad in adapter_dict:
            if ad in seq:
                ad_count[ad] += 1
                # ad_pos[ad].add(seq.index(ad))
                pos = seq.index(ad)
                ad_pos[ad][pos] += 1
                # if pos < 100:
                #     print(f"len: {len(seq)}  pos: {pos}")
                #     print(seq)
    for ad in ad_count:
        print(f"{adapter_dict[ad]} ({ad}) is present {ad_count[ad]} times in read-fq")
        # print(ad_pos[ad])
        # print(ad_pos)
        # ad_pos_sort = sorted(ad_pos.items())
        for n in range(max(ad_pos[ad])+1):
            col = int((ad_pos[ad][n]/20))
            if col > 0:
                print(f"{n}\t{'*' * col}")
    if extracted:
        os.remove(extracted)
    return ad_count


if __name__ == "__main__":
    adapter_fa = sys.argv[1]
    read_fa = sys.argv[2]
    start = perf_counter()
    ad_count = chk_ilmn_reads (read_fa, adapter_fa)
    stop = perf_counter()
    print(f"Runtime: {stop-start}")