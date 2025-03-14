import json

configfile: "config/NanoporeOnlyVirus.yaml"

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
        "python scripts/preprocessing.py {input} {wildcards.experiment} {output} 2>&1 > {log}"


rule get_filter_params:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    params:
        # target_cov="80"
        target_cov=config["TargetCoverage"],
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
        fastq="results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq",
        params_json="results/{experiment}/{barcode}/filt_params.json"
    params:
        # minlen="2000"  # It would be best if these params would be automatically determined with select_filter_params_by_cov_sm.py or from config
        minlen= lambda wildcards, input: json.load(open(input.params_json))["len"]
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sizefilt.fastq"
    log:
        "logs/{experiment}/{barcode}/size_filter.log"
    conda:
        "envs/bbmap.yaml"
    shell:
        "bbduk.sh in={input.fastq} out={output} minlen={params.minlen} 2>&1 > {log}"


rule filter_readqual:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sizefilt.fastq",
        params_json="results/{experiment}/{barcode}/filt_params.json"
    params:
        # minqual="10"  # It would be best if these params would be automatically determined with select_filter_params_by_cov_sm.py or from config
        minqual= lambda wildcards, input: json.load(open(input.params_json))["qual"]
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq"
    log:
        "logs/{experiment}/{barcode}/qual_filter.log"
    conda:
        "envs/bbmap.yaml"
    shell:
        "bbduk.sh in={input} out={output} maq={params.minqual} 1> {log}"


rule porechop_barcodes_and_adapters:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
    log:
        "logs/{experiment}/{barcode}/porechop.log"
    conda:
        "envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} 2>&1 > {log}"


rule assemble_canu:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/canu"),
        contigs="results/{experiment}/{barcode}/canu/{experiment}_{barcode}_canu.contigs.fasta"
    params:
        # genomeSize=230000,
        genomeSize=config["GenomeSize"],
        useGrid="False"
    # threads:
    #     16
    log:
        "logs/{experiment}/{barcode}/canu.log"
    conda:
        "envs/canu.yaml"
    shell:
        "canu -nanopore {input} -p {wildcards.experiment}_{wildcards.barcode}_canu -d {output.outdir} genomeSize={params.genomeSize} useGrid={params.useGrid} 2>&1 > {log}"


# rule polish_canu_nextpolish:
#     input:
#         reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta",
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
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta",
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
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/flye"),
        outfile="results/{experiment}/{barcode}/flye/assembly.fasta"
    params:
        # genomeSize=230000,
        genomeSize=config["GenomeSize"],
        useGrid="False"
    threads:
        16
    log:
        "logs/{experiment}/{barcode}/canu.log"
    conda:
        "envs/flye.yaml"
    shell:
        "flye --nano-hq {input} -o {output.outdir} -g {params.genomeSize} -t {threads} 2>&1 > {log}"


rule polish_flye_medaka:
    input:
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta",
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

    