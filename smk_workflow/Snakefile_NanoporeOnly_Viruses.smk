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


rule filter_out_host_reads: # design in such a way, that if the assembly is from a pure culture an empty mock genome can be used, so nothing will be filtered
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq"
    params:
        host_gen_ref=config["HostRefGenome"]
    log:
        "logs/{experiment}/{barcode}/filter_out_host_reads.log"
    conda:
        "envs/minimap2.yaml"
    shell:
        "python scripts/filter_out_host_reads.py -q {input} -hgr {params.host_gen_ref} -o {output} 2>&1 > {log}"


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
        minqual= lambda wildcards, input: json.load(open(input.params_json))["qual"]
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq"
    log:
        "logs/{experiment}/{barcode}/qual_filter.log"
    conda:
        "envs/bbmap.yaml"
    shell:
        "bbduk.sh in={input} out={output} maq={params.minqual} 1> {log}"


rule get_longest_reads:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq" # use non-host
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_nonhost.fastq"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq.longestx"
    params:
        longest_cov=config["LongestReadsCov"],
        genomeSize=config["GenomeSize"]
    log:
        "logs/{experiment}/{barcode}/get_longest_reads.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/get_x_cov_longest_reads.py -cov {params.longest_cov} -g {params.genomeSize} -r {input} 2>&1 > {log}"


rule inject_longest_reads_into_filtered:
    input:
        sqfilt="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq",
        longest="results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq.longestx"
    output:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fastq"
    log:
        "logs/{experiment}/{barcode}/inject_longest_reads_into_filtered.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/inject_longest_reads_into_filt.py -a {input.sqfilt} -i {input.longest} 2>&1 > {log}"
    
    

rule porechop_barcodes_and_adapters:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fastq"
    output:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fasta"
    log:
        "logs/{experiment}/{barcode}/porechop.log"
    conda:
        "envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} 2>&1 > {log}"


rule convert_multiline_fasta:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fasta"
    output:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        
    log:
        "logs/{experiment}/{barcode}/porechop.log"
    conda:
        "envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} 2>&1 > {log}"


rule assemble_canu:
    input:
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/canu"),
        contigs="results/{experiment}/{barcode}/canu/{experiment}_{barcode}_canu.contigs.fasta"
    params:
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
        # "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/flye"),
        outfile="results/{experiment}/{barcode}/flye/assembly.fasta"
    params:
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
        # reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta",
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.pluslong.fasta",
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