import os
import json
import glob


configfile: "config/Snakefile_NanoporeOnly_Metagenomic_Viruses_Bacteria.yaml"

EXPERIMENT = config["Experiment"]
BARCODES = config["Barcodes"]
ASSEMBLER = config["Assembler"]


rule all:
    input:
        expand("results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}/consensus.oriented.fasta", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLER),
        # expand("results/{experiment}/{barcode}/busco/medaka_{assembler}/logs/busco.log", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLER),
        expand("results/{experiment}/{barcode}/minikraken_{assembler}/custom_summary.tsv", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLER),
        expand("results/{experiment}/{barcode}/minikraken_reads/classifications_nonhost_reads.report", experiment=EXPERIMENT, barcode=BARCODES)

# concatentation should not be done in this step (only decompressing) as it is
# now, because the concatenation is done after filtering of host-reads now, which
# is (currently) done on the single fastq level, not with the catenated big file
# due to I/O slow down with huge files
rule preprocessing:
    input:
        "data/{experiment}/{barcode}"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    log:
        "logs/{experiment}/{barcode}/preprocessing.log"
    conda:
        "envs/python3.yaml"
    shell:
        # think about modifiyng pre-proc script to contain argparse and options -unzip and -cat, because concatenation may not be needed anymore because of individual filtering 
        "python scripts/preprocessing.py {input} {wildcards.experiment} {output} 2>&1 > {log}"


# def get_input_fastq(wildcards):
#     # might be amended to do the same for fasta to be more flexible
#     read_dir = f"data/{wildcards.experiment}/{wildcards.barcode}"
#     return glob.glob(os.path.join(read_dir, "*.fastq"))

# checkpoint filter_out_host_reads:
#     output:
#         directory("results/{experiment}/{barcode}/non_host_reads")
#     shell:
        

rule filter_out_host_reads: # design in such a way, that if the assembly is from a pure culture an empty mock genome can be used, so nothing will be filtered
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
        # get_input_fastq
        # "results/{experiment}/{barcode}/{fastq}.fastq"
        "data/{experiment}/{barcode}"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq",
        # individual = "results/{experiment}/{barcode}/{fastq}_nonhost.fastq"
    params:
        host_gen_ref=config["HostRefGenome"],
        outdir = "results/{experiment}/{barcode}/nonhost_fastq"
    threads:
        # it would be best to implement a dynamic solution in the script which detects 
        # the amount of RAM and decides based on that how many workers to initiate
        # workflow.cores
        12
    log:
        "logs/{experiment}/{barcode}/filter_out_host_reads.log"
    conda:
        "envs/minimap2.yaml"
    shell:
        "python scripts/filter_out_host_reads_multi.py --cores {threads} -d {input} -hgr {params.host_gen_ref} -o {params.outdir} 2>&1 > {log}"


rule kraken2_classification_reads:
    input:
        fastq = "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq"
    output:
        raw = "results/{experiment}/{barcode}/minikraken_reads/classifications_nonhost_reads.kraken",
        report = "results/{experiment}/{barcode}/minikraken_reads/classifications_nonhost_reads.report"
    params:
        db = config["KrakenDB"]
    log:
        "logs/{experiment}/{barcode}/kraken2/minikraken_reads/kraken2.log"
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {params.db} \
                --output {output.raw} \
                --report {output.report} \
                {input.fastq}
        """


rule get_filter_params:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq"
    params:
        target_cov=config["TargetCoverageAsm"],
        genomeSize=config["GenomeSize"]
    output:
        "results/{experiment}/{barcode}/filt_params.json"
    log:
        "logs/{experiment}/{barcode}/get_filter_params.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/select_filter_params_by_cov_sm.py -c {params.target_cov} -s {params.genomeSize} --allow_reduction_all -f {input} 2>&1 > {log}"
    

rule filter_readlength:
    input:
        # fastq="results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq",
        fastq="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq",
        params_json="results/{experiment}/{barcode}/filt_params.json"
    params:
        minlen= lambda wildcards, input: json.load(open(input.params_json))["len"]
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sizefilt.fastq"
    log:
        "logs/{experiment}/{barcode}/size_filter.log"
    conda:
        "envs/bbmap.yaml"
    shell:
        "bbduk.sh in={input.fastq} out={output} minlen={params.minlen} 2>&1 > {log}"


rule filter_readqual:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sizefilt.fastq",
        params_json="results/{experiment}/{barcode}/filt_params.json"
    params:
        minqual= lambda wildcards, input: json.load(open(input.params_json))["qual"]
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.fastq"
    log:
        "logs/{experiment}/{barcode}/qual_filter.log"
    conda:
        "envs/bbmap.yaml"
    shell:
        "bbduk.sh in={input} out={output} maq={params.minqual} 1> {log}"


rule get_longest_reads:
    input:
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq",
        # the below is added here only to prevent this rule from running in parallel with 
        # get_filter_params, because both scripts access the same file and use next()
        # which causes errors/race conditions in both processes
        start_flag="results/{experiment}/{barcode}/filt_params.json" 
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq.longestx"
    params:
        longest_cov=config["LongestReadsCov"],
        genomeSize=config["GenomeSize"]
    log:
        "logs/{experiment}/{barcode}/get_longest_reads.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/get_x_cov_longest_reads.py -cov {params.longest_cov} -g {params.genomeSize} -r {input.reads} 2>&1 > {log}"


rule inject_longest_reads_into_filtered:
    input:
        sqfilt="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.fastq",
        longest="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq.longestx"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fastq"
    log:
        "logs/{experiment}/{barcode}/inject_longest_reads_into_filtered.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/inject_longest_reads_into_filt.py -a {input.sqfilt} -i {input.longest} 2>&1 > {log}"
    

rule porechop_barcodes_and_adapters:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fastq"
    output:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fasta"
    log:
        "logs/{experiment}/{barcode}/porechop.log"
    conda:
        "envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} 2>&1 > {log}"


# rule convert_multiline_fasta:
#     input:
#         "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fasta"
#     output:
#         # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.fasta"   



rule assemble_canu:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/canu"),
        contigs="results/{experiment}/{barcode}/canu/{experiment}_{barcode}_canu.contigs.fasta"
    params:
        genomeSize=config["GenomeSize"],
        useGrid="False",
        minInputCoverage=config["CanuMinInput"],
        stopOnLowCoverage=config["CanuStopOnCov"]

    # threads:
    #     16
    log:
        "logs/{experiment}/{barcode}/canu.log"
    conda:
        "envs/canu.yaml"
    shell:
        "canu -nanopore {input} -p {wildcards.experiment}_{wildcards.barcode}_canu -d {output.outdir} genomeSize={params.genomeSize} useGrid={params.useGrid} minInputCoverage={params.minInputCoverage} stopOnLowCoverage={params.stopOnLowCoverage} 2>&1 > {log}"


# rule polish_canu_nextpolish:
#     input:
#         reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.fasta",
#         draft_asm="results/{experiment}/{barcode}/canu/{experiment}_{barcode}_canu.contigs.fasta"
#     output:
#         outdir=directory("results/{experiment}/{barcode}/medaka_canu"),
#         consensus="results/{experiment}/{barcode}/medaka_canu/consensus.fasta"
#     log:
#         "logs/{experiment}/{barcode}/medaka_canu.log"
#     threads:
#         16
#     conda:
#         "envs/medaka.yaml"
#     shell:
#         "medaka_consensus -i {input.reads} -d {input.draft_asm} -o {output.outdir} -t {threads} 2>&1 > {log}"


rule polish_canu_medaka:
    input:
        # reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.fasta",
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fasta",
        draft_asm="results/{experiment}/{barcode}/canu/{experiment}_{barcode}_canu.contigs.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/medaka_canu"),
        consensus="results/{experiment}/{barcode}/medaka_canu/consensus.fasta"
    log:
        "logs/{experiment}/{barcode}/medaka_canu.log"
    threads:
        16
    conda:
        "envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.reads} -d {input.draft_asm} -o {output.outdir} -t {threads} 2>&1 > {log}"

rule assemble_flye:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/flye"),
        outfile="results/{experiment}/{barcode}/flye/assembly.fasta"
    params:
        genomeSize=config["GenomeSize"],
        useGrid="False"
    threads:
        16
    log:
        "logs/{experiment}/{barcode}/flye.log"
    conda:
        "envs/flye.yaml"
    shell:
        "flye --nano-hq {input} -o {output.outdir} -g {params.genomeSize} -t {threads} 2>&1 > {log}" # consider using --meta 

# consider adding ne rule for metaFlye assembler (separate tool derived from Flye)

rule polish_flye_medaka:
    input:
        # reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta",
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost_sqfilt.pluslong.fasta",
        draft_asm="results/{experiment}/{barcode}/flye/assembly.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/medaka_flye"),
        consensus="results/{experiment}/{barcode}/medaka_flye/consensus.fasta"
    log:
        "logs/{experiment}/{barcode}/medaka_flye.log"
    threads:
        16
    conda:
        "envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.reads} -d {input.draft_asm} -o {output.outdir} -t {threads} 2>&1 > {log}"


#  rule generate_medaka_report:
#      input:
#          "results/{experiment}/{barcode}/medaka_flye/consensus.fasta"
#      output:

rule circlator_fixstart:
    input:
        "results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta"
    output:
        # "results/{experiment}/{barcode}/circl_fixstart/{assembly}.oriented.fasta"
        "results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}/consensus.oriented.fasta"
    params:
        out_prefix = "results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}/consensus.oriented"
    conda:
        "envs/circlator.yaml"
    log:
        "logs/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}/circlator_fixstart.log"
    shell:
        "circlator fixstart {input} {params.out_prefix} 2>&1 > {log}"

# create a rule that splits the assemblies into separate contigs, so busco can 
# assess them individually using automatic lineage detection; either split
# them blindly by headers or try a tool like kraken that tries to identify their
# lineage/species identity first

rule kraken2_classification_assembly:
    input:
        fasta = "results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta"
    output:
        raw = "results/{experiment}/{barcode}/minikraken_{assembler}/classifications.kraken",
        report = "results/{experiment}/{barcode}/minikraken_{assembler}/classifications.report"
    params:
        db = config["KrakenDB"]
    log:
        "logs/{experiment}/{barcode}/kraken2/minikraken_{assembler}/kraken2.log"
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {params.db} \
                --output {output.raw} \
                --report {output.report} \
                {input.fasta}
        """
# snakemake --cores 32 --use-conda results/SM0037_CMV_Voigt/barcode01/minikraken_flye/classifications.kraken


rule kraken_assembly_post_processing:
    input:
        "results/{experiment}/{barcode}/minikraken_{assembler}/classifications.report"
    output:
        "results/{experiment}/{barcode}/minikraken_{assembler}/custom_summary.tsv"
    params:
        kraken_outdir="results/{experiment}/{barcode}/minikraken_{assembler}"
    log:
        "logs/{experiment}/{barcode}/minikraken_{assembler}_post_processing/kraken_post_processing.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/kraken_assembly_classification_post_processing.py {params.kraken_outdir} 2>&1 > {log}"


rule busco:
    input:
        "results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta"
    output:
        # "results/{experiment}/{barcode}/medaka_{assembler}/"
        "results/{experiment}/{barcode}/busco/medaka_{assembler}/logs/busco.log"
    params:
        # out="results/{experiment}/{barcode}/medaka_{assembler}/busco_out"
        out="results/{experiment}/{barcode}/busco/medaka_{assembler}"
    conda:
        "envs/busco.yaml"
    log:
        "logs/{experiment}/{barcode}/busco/medaka_{assembler}/busco.log"
    shell: # use batch mode for split assembly
         # -m is the mode flag => run in genome mode (altern.: "transcriptome" and "proteins")
        "busco -i {input} -o {params.out} -m genome --auto-lineage --cpu 16 2>&1 > {log}"
        # "busco -i metagenome.fasta -m genome -l bacteria_odb12,viruses_odb12 --cpu 16"
    
# maybe add rule for creation of busco reports either like in big bacterial workflow or using buscomp tool for compilation of busco results
# maybe add rule and script that exclude contigs below a set min-length, to filter out truncated and/or badly assmebled contigs