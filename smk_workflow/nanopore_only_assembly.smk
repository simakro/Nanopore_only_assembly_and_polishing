

rule preprocessing:
    input:
        "data/{experiment}/{barcode}"
    output:
        "data/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    log:
        "logs/{experiment}/{barcode}/preprocessing.log"
    conda:
        "../envs/python3.yaml"
    shell:
        "python /scripts/preprocessing.py {input} {wildcards.experiment}"

        



# rule filter_readlength:
#     input:



# def filter_readlength(cat_fq, min_len):
#     infile = f"in={cat_fq}"
#     path, fname = os.path.split(cat_fq)
#     outname= os.path.join(path, f"{'.'.join(fname.split('.')[:-1])}_{int(min_len/1000)}kb.fastq")
#     outfile = f"out={outname}"
#     cutoff = f"minlen={min_len}"
#     sp.run(["bbduk.sh", infile, outfile, cutoff])
#     return outname


# def filter_quality(len_filtered, min_qual):
#     infile = f"in={len_filtered}"
#     path, fname = os.path.split(len_filtered)
#     outname = os.path.join(path, f"{'.'.join(fname.split('.')[:-1])}_Q{int(min_qual)}.fastq")
#     outfile = f"out={outname}"
#     cutoff = f"maq={min_qual}"
#     sp.run(["bbduk.sh", infile, outfile, cutoff])
#     return outname


# def porechop_barcodes_and_adapters(ql_filt_reads):
#     outname = ".".join(ql_filt_reads.split(".")[:-1]) + ".fasta"
#     sp.run([
#         "porechop",
#         "-i", ql_filt_reads,
#         "-o", outname
#         ])
#     return outname


# # def guppy_barcoder_trimming(ql_filt_reads):
# #     root, leaf = os.path.split(ql_filt_reads)
# #     savepath = os.path.join(root, "bc_trimmed")
# #     sp.run([
# #         "guppy_barcoder",
# #         "--barcode_kits", "EXP-NBD104",
# #         "--enable_trim_barcodes",
# #         "--trim_adapters",
# #         "--save_path", savepath,
# #         "-i", ql_filt_reads,
# #         ])
# #     return outname

# def assemble_canu(preprocessed_reads):
#     root, leaf = os.path.split(preprocessed_reads)
#     asm_name = ".".join(leaf.split(".")[:1]) + "_assembly_canu"
#     asm_dir = os.path.join(root, asm_name)
#     prefix = "_".join(leaf.split("_")[:2]) + "_"
#     sp.run([
#         "canu",
#         f"-nanopore", preprocessed_reads,
#         "-p", prefix,
#         "-d", asm_dir,
#         "genomeSize=230000",
#         "useGrid=False"
#         ])


# def polish_pepper():
#     pass

# def polish_medaka():
#     pass

# def polish_homopolish()
#     pass

# # def polish_nanopolish():
# #     pass



# if __name__ == "__main__":
#     read_folder = sys.argv[1]
#     experiment = sys.argv[2]
#     min_readlen = int(sys.argv[3])
#     min_qual = int(sys.argv[4])
#     decompress_all(read_folder)
#     cat_fq = concatenate_fastq(read_folder, experiment)
#     len_filt = filter_readlength(cat_fq, min_readlen)
#     q_l_flt = filter_quality(len_filt, min_qual)
#     pp_reads = porechop_barcodes_and_adapters(q_l_flt)
#     # pp_reads = guppy_barcoder_trimming(q_l_flt)
#     assemble_canu(pp_reads)

