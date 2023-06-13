import os
import sys
import gzip


def decompress_all(read_folder):
    files = os.scandir(read_folder)
    for f in files:
        if f.path.endswith(".gz"):
            outfile = ".".join(f.path.split(".")[:-1])
            with open(outfile, "w") as out, gzip.open(f.path, "r") as data:
                for line in data:
                    out.write(line.decode("utf-8"))
            os.remove(f.path)


def concatenate_fastq(read_folder, exp_prefix):
    files = list(os.scandir(read_folder))
    outname = f"{exp_prefix}_{os.path.split(read_folder)[-1]}_all.fastq"
    outfile = os.path.join(read_folder, outname)
    if os.path.exists(outfile):
        print("Concatenated file already exists")
    else:
        with open(outfile, "w") as out:
            for f in files:
                print(f.path)
                if f.path.endswith(".fastq") and not f.path.endswith("_all.fastq"):
                    with open(f.path, "r") as fq:
                        data = fq.read()
                        for line in data:
                            out.write(line)
                    os.remove(f.path)
    return outfile

if __name__ == "__main__":
    read_folder = sys.argv[1]
    exp_prefix = sys.argv[2]
    decompress_all(read_folder)
    concatenate_fastq(read_folder, exp_prefix)