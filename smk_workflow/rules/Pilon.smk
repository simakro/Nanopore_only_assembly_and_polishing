
include: "Common.smk"


rule bwa_ilmn_to_rawasm:
    input:
        raw_asm="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
        ilmn_reads=get_clipped_ilmn_reads
    output:
        bam_sort="results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam.sort",
        sam_file=temp("idx1_aln_ilm_{experiment}_{barcode}.sam"),
        bam=temp("results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam")
    params:
        prefix="idx1_medaka_{assembler}_{experiment}_{barcode}",
        outdir="results/{experiment}/{barcode}/medaka_{assembler}_bwa"#,
        # bam="results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam"
        # sam_file="idx1_aln_ilm_{experiment}_{barcode}.sam",
    threads:
        16
    conda:
        "../envs/bwa-mem2.yaml"
    log:
        "logs/{experiment}/{barcode}/bwa-mem-1_medaka_{assembler}.log"
    shell:
        "mkdir -p {params.outdir} && "
        "bwa-mem2 index {input.raw_asm} -p {params.prefix} && "
        "bwa-mem2 mem -t {threads} {params.prefix} {input.ilmn_reads} > {output.sam_file} && "
        "samtools view -hbS {output.sam_file} > {output.bam} && "
        "samtools sort {output.bam} > {output.bam_sort} && "
        "samtools index {output.bam_sort} && "
        "mv {params.prefix}* results/{wildcards.experiment}/{wildcards.barcode}/medaka_{wildcards.assembler} && "
        "mv {output.sam_file} {params.outdir} 2>&1 > {log}"


rule pilon_raw_asm:
    input:
        genome="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
        bam="results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam.sort"
    output:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon1/medaka_{assembler}_pilon1.fasta",
        flag="results/{experiment}/{barcode}/medaka_{assembler}_pilon1/medaka_{assembler}_pilon1.complete.flag"
    params:
        prefix="medaka_{assembler}_pilon1",
        outdir="results/{experiment}/{barcode}/medaka_{assembler}_pilon1"
    conda:
        "../envs/pilon.yaml"
    log:
        "logs/{experiment}/{barcode}/pilon-1_medaka_{assembler}.log"
    shell:
        "pilon -Xmx16G --genome {input.genome} --bam {input.bam} --output {params.prefix} --outdir {params.outdir} 2>&1 > {log} && "
        "touch {output.flag}"


rule pilon_iteration_2:
    input:
        asm = "results/{experiment}/{barcode}/medaka_{assembler}_pilon1/medaka_{assembler}_pilon1.fasta",
        ilmn_reads=get_clipped_ilmn_reads
    output:
        sam=temp("results/{experiment}/{barcode}/medaka_{assembler}_pilon2/medaka_{assembler}_ilmn2pilon1.sam"),
        bam=temp("results/{experiment}/{barcode}/medaka_{assembler}_pilon2/idx2_medaka_{assembler}_ilmn2pilon1.bam"),
        sort="results/{experiment}/{barcode}/medaka_{assembler}_pilon2/idx2_medaka_{assembler}_ilmn2pilon1.bam.sort",
        asm="results/{experiment}/{barcode}/medaka_{assembler}_pilon2/medaka_{assembler}_pilon2.fasta"
    params:
        # threads=8,
        outdir="results/{experiment}/{barcode}/medaka_{assembler}_pilon2",
        bwa_prefix="results/{experiment}/{barcode}/medaka_{assembler}_pilon2/idx2_medaka_{assembler}_pilon1",
        pilon_prefix="medaka_{assembler}_pilon2",
        # sam="results/{experiment}/{barcode}/medaka_{assembler}_pilon2/medaka_{assembler}_ilmn2pilon1.sam",
        # bam="results/{experiment}/{barcode}/medaka_{assembler}_pilon2/idx2_medaka_{assembler}_ilmn2pilon1.bam"
    threads:
        8
    conda:
        "../envs/pilon_iteration.yaml"
    log:
        "logs/{experiment}/{barcode}/medaka_{assembler}_pilon2.log"
    shell:
        "mkdir -p {params.outdir} && "
        "bwa-mem2 index {input.asm} -p {params.bwa_prefix} && "
        "bwa-mem2 mem -t {params.threads} {params.bwa_prefix} {input.ilmn_reads} > {output.sam}  && "
        "samtools view -hbS {output.sam} > {output.bam} && "
        "samtools sort {output.bam} > {output.sort} && "
        "samtools index {output.sort} && "
        "pilon -Xmx16G --genome {input.asm} --bam {output.sort} --output {params.pilon_prefix} --outdir {params.outdir} 2>&1 > {log}"
        