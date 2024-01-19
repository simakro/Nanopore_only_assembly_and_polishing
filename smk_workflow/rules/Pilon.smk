
include: "Common.smk"


rule bwa_ilmn_to_rawasm:
    input:
        raw_asm="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
        # ilmn_r1,ilmn_r2=get_ilmn_reads
        # flag="results/{experiment}/{barcode}/ilmn_fastqc_after/fastqc_complete.flag",
        ilmn_reads=get_clipped_ilmn_reads
    output:
        "results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam.sort"#,
        # outdir=directory("results/{experiment}/{barcode}/medaka_{assembler}/bwa")
    params:
        threads=16,
        prefix="idx1_medaka_{assembler}_{experiment}_{barcode}",
        sam_file="idx1_aln_ilm_{experiment}_{barcode}.sam",
        outdir="results/{experiment}/{barcode}/medaka_{assembler}_bwa",
        bam="results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam"
    conda:
        "../envs/bwa-mem2.yaml"
    log:
        "logs/{experiment}/{barcode}/bwa-mem-1_medaka_{assembler}.log"
    shell:
        "mkdir -p {params.outdir} && "
        "bwa-mem2 index {input.raw_asm} -p {params.prefix} && "
        "bwa-mem2 mem -t {params.threads} {params.prefix} {input.ilmn_reads} > {params.sam_file} && " # {input.ilmn_r1} {input.ilmn_r2}
        "samtools view -hbS {params.sam_file} > {params.bam} && "
        "samtools sort {params.bam} > {output} && "
        "samtools index {output} && "
        # "mv {params.prefix}.bam results/{wildcards.experiment}/{wildcards.barcode}/medaka_{wildcards.assembler} && " # this is just to test if the pilon WARNING that bam.bai index is older then BAM is related to moving the files around and getting new time stamps
        "mv {params.prefix}* results/{wildcards.experiment}/{wildcards.barcode}/medaka_{wildcards.assembler} && "
        "mv {params.sam_file} {params.outdir} 2>&1 > {log}"


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



# rule bwa_mem_iteration:
#     input:
#         # raw_asm="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
#         # "results/{experiment}/{barcode}/medaka_{assembler}_pilon{iteration}/medaka_{assembler}_pilon{iteration}.fasta"
#         # genome=get_genome_prev_iter,
#         # genome = lambda wildcards: "results/{experiment}/{barcode}/medaka_{assembler}_pilon{0}/medaka_{assembler}_pilon{0}.fasta".format(int(wildcards.iteration)-1),
#         genome = lambda wildcards: "results/{0}/{1}/medaka_{2}_pilon{3}/medaka_{2}_pilon{3}.fasta".format(wildcards.experiment, wildcards.barcode, wildcards.assembler, int(wildcards.iteration-1)),
#         # ilmn_r1,ilmn_r2=get_ilmn_reads
#         # ilmn_reads=get_clipped_ilmn_reads
#     output:
#         "results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx{iteration}_aln_ilm_{experiment}_{barcode}.bam.sort"#,
#         # outdir=directory("results/{experiment}/{barcode}/medaka_{assembler}/bwa")
#     params:
#         ilmn_reads=get_clipped_ilmn_reads,
#         threads=16,
#         prefix="idx{iteration}_medaka_{assembler}_{experiment}_{barcode}",
#         sam_file="idx{iteration}_aln_ilm_{experiment}_{barcode}.sam",
#         outdir="results/{experiment}/{barcode}/medaka_{assembler}_bwa",
#         bam="results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx{iteration}_aln_ilm_{experiment}_{barcode}.bam"
#     conda:
#         "../envs/bwa-mem2.yaml"
#     log:
#         "logs/{experiment}/{barcode}/bwa-mem-{iteration}_medaka_{assembler}.log"
#     shell:
#         "mkdir -p {params.outdir} && "
#         "bwa-mem2 index {input.raw_asm} -p {params.prefix} && "
#         "bwa-mem2 mem -t {params.threads} {params.prefix} {params.ilmn_reads} > {params.sam_file} && " # {input.ilmn_r1} {input.ilmn_r2}
#         "samtools view -hbS {params.sam_file} > {params.bam} && "
#         "samtools sort {params.bam} > {output} && "
#         "samtools index {output} && "
#         # "mv {params.prefix}.bam results/{wildcards.experiment}/{wildcards.barcode}/medaka_{wildcards.assembler} && " # this is just to test if the pilon WARNING that bam.bai index is older then BAM is related to moving the files around and getting new time stamps
#         "mv {params.prefix}* results/{wildcards.experiment}/{wildcards.barcode}/medaka_{wildcards.assembler} && "
#         "mv {params.sam_file} {params.outdir} 2>&1 > {log}"


# rule pilon_iteration:
#     input:
#         # genome="results/{experiment}/{barcode}/medaka_{assembler}_pilon/medaka_{assembler}_pilon1.fasta"
#         # genome=get_genome_prev_iter,
#         genome = lambda wildcards: "results/{0}/{1}/medaka_{2}_pilon{3}/medaka_{2}_pilon{3}.fasta".format(wildcards.experiment, wildcards.barcode, wildcards.assembler, int(wildcards.iteration-1)),
#         # genome="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
#         bam="results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx{iteration}_aln_ilm_{experiment}_{barcode}.bam.sort"
#     output:
#         "results/{experiment}/{barcode}/medaka_{assembler}_pilon{iteration}/medaka_{assembler}_pilon{iteration}.fasta"
#     params:
#         prefix="medaka_{assembler}_pilon{iteration}",
#         outdir="results/{experiment}/{barcode}/medaka_{assembler}_pilon{iteration}"
#     conda:
#         "../envs/pilon.yaml"
#     log:
#         "logs/{experiment}/{barcode}/pilon-{iteration}_medaka_{assembler}.log"
#     shell:
#         "pilon -Xmx16G --genome {input.genome} --bam {input.bam} --output {params.prefix} --outdir {params.outdir} 2>&1 > {log}"