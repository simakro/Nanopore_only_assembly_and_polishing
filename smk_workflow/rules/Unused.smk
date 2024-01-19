# potential tools for adapter clipping/trimming: Scythe, Trimmomatic, cutadapt
# adapter clipping tested and discarded: trim_galore (wrapping cutadapt together with fastqc)

# potential tools for quality trimming: Sickle, Trimmomatic, 

# rule clip_adapters_scythe:
#     input:
#         "results/{experiment}/{barcode}/ilmn_fastqc/fastqc_complete.flag", # this is not required as input for fastqc, but supplies the necessary wildcards and ensures 
#         ilmn_reads=get_ilmn_reads
#     output:
#     params:
#         outdir="results/{experiment}/{barcode}/trim_scythe",
#         adapters=""
#     conda:
#         "envs/scythe.yaml"
#     shell:
#         "scythe -a {params.adapters} -o {input.ilmn_reads}"
#         # "fastqc {input.ilmn_reads} -o {params.outdir}"


# rule homopolish: # Reference guided/based polishing without read information and may conceal variants! Assumes close homology (hence the name) between ref and sample.
#     input:
#         asm="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
#         # ref="resources/KP202868.1_Murid_herpesvirus_8_isolate_Berlin.fa"
#         ref= "resources/{experiment}/{reference}"
#     output:
#         outdir = directory("results/{experiment}/{barcode}/homopolish_medaka_{assembler}"),
#         outfile= "results/{experiment}/{barcode}/homopolish_medaka_{assembler}/consensus_homopolished.fasta"
#     log:
#         "logs/{experiment}/{barcode}/medaka_{assembler}_homopolish.log"
#     threads:
#         16
#     conda:
#         "envs/homopolish.yaml"
#     shell:
#         "homopolish polish -a {input.asm} -l {input.ref} -m R9.4.pkl -o {output.outdir} -t {threads} 2>&1 > {log}"