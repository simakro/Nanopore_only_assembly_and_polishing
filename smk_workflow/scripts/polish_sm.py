import argparse
from time import perf_counter


def get_args():
    parser = argparse.ArgumentParser()
    
    required_args = parser.add_argument_group("Required Arguments")
    required_args.add_argument(
        "-a", "--asm", required=True,
        help="Path to (long read) genome assembly"
    )
    required_args.add_argument(
        "-r", "--short_reads", required=True,
        help="Path to high quality (short) reads"
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


# def load_adapter(adapter_fa):
#     adapter_dict = dict()
#     with open(adapter_fa, "r") as fa:
#         for line in fa:
#             if line.startswith(">"):
#                 header = line.strip()
#                 seq = next(fa).strip().upper()
#                 adapter_dict[seq] = header
#     return adapter_dict


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


def sieve_mismatching(args):
    asm = open(args.asm,"r").read()
    read_fq = args.short_reads
    if read_fq.endswith(".gz"):
        print("decompressing")
        read_fq = extract_gz(read_fq, keep_gz=True)
        extracted = read_fq
    fq_gen = sl_fastq_generator(read_fq)
    all_count = 0
    mismatch_count = 0
    with open(read_fq + ".mismatches", "w") as mm:
        for rid,seq,plus,qual in fq_gen:
            all_count += 1
            if seq not in asm:
                mismatch_count += 1
                read_data = [rid,seq,plus,qual]
                for item in read_data:
                    mm.write(f"{item}\n")
    print(f"Analyzed {all_count} short reads.")
    perc = (mismatch_count/all_count)*100
    print(f"Found {mismatch_count} mismatching short reads ({per}%).")


if __name__ == "__main__":
    start = perf_counter()
    args = get_args()
    sieve_mismatching(args)
    end = perf_counter()
    print(f"Runtime: {end-start} sec")




