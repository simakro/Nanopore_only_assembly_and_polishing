import os
import sys
import gzip


def ensure_naming_convention(read_folder):
    files = list(os.scandir(read_folder))
    expected_ext = ["fastq"]
    unexpected_files = []
    for f in files:
        fnm = f.path
        ext = fnm.split(".")[-1]
        if ext not in expected_ext:
            unexpected_files.append(fnm)
        else:
            if f.path.endswith("_all.fastq"):
                print('Input fastq files may not end on "_all.fastq"')
                print(f.path)
                print("violates naming convention")
                sys.exit()
    if not len(unexpected_files):
        pass
    else:
        print(
            f"The following unexpected files were found in input folder {unexpected_files}"
        )


def decompress_all(read_folder):
    files = os.scandir(read_folder)
    for f in files:
        if f.path.endswith(".gz"):
            outfile = ".".join(f.path.split(".")[:-1])
            with open(outfile, "w") as out, gzip.open(f.path, "r") as data:
                for line in data:
                    out.write(line.decode("utf-8"))
            os.remove(f.path)


def concatenate_fastq(read_folder, exp_prefix, output):
    files = list(os.scandir(read_folder))
    # outname = f"{exp_prefix}_{os.path.split(read_folder)[-1]}_all.fastq"
    # outfile = os.path.join(read_folder, outname)
    if os.path.exists(output):
        print("Concatenated file already exists")
    else:
        with open(output, "w") as out:
            for f in files:
                print(f.path)
                if f.path.endswith(".fastq") and not f.path.endswith("_all.fastq"):
                    with open(f.path, "r") as fq:
                        data = fq.read()
                        for line in data:
                            out.write(line)
    return output

if __name__ == "__main__":
    read_folder = sys.argv[1]
    exp_prefix = sys.argv[2]
    output = sys.argv[3]
    decompress_all(read_folder)
    ensure_naming_convention(read_folder)
    concatenate_fastq(read_folder, exp_prefix, output)
