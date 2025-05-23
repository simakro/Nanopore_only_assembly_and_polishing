import os
import glob
import json
from collections import defaultdict


include: "rules/Pilon.smk"
include: "rules/Common.smk"


EXPERIMENT = "SM0026_bacterial_genomes_KP"
BARCODES = [
    "barcode01",
    "barcode02",
    "barcode03",
    "barcode04",
    "barcode05",
    "barcode06",
    "barcode07",
    "barcode08",
    "barcode09",
    "barcode10",
    "barcode11",
    "barcode12",
    "barcode13",
    ]

with open(os.path.join("data", EXPERIMENT, "Sample_Info.json")) as sinfo:
    SAMPLE_INFO = {k:v for k,v in json.load(sinfo).items() if k in BARCODES}
    # print("SAMPLE_INFO", SAMPLE_INFO)

with open(os.path.join("resources", "GTDBTK_classification_info.json")) as gcinfo:
    GTDB_TRANSL = json.load(gcinfo)

ASSEMBLERS = ["flye"] # , "canu"
GENOME_SIZE = 3000000
COVERAGE = 50
REFERENCES = [".".join(SAMPLE_INFO[barcode]["ref_fallback"].split(".")[:-1]) for barcode in SAMPLE_INFO]
BUSCO_GENUS = "bacteria"
PILON_ITERATIONS = 1

rule all:
    input:
        expand("results/{experiment}/{barcode}/prokka_medaka_{assembler}_protgbk_class/{experiment}_{barcode}_{assembler}_prokka_protgbk_class.tsv", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/busco_medaka_{assembler}_pilon3/run_{genus}_odb10/short_summary.specific.{barcode}.txt", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS, genus=BUSCO_GENUS),
        expand("results/{experiment}/busco_graph/{assembler}/busco_figure.png", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_bwa/idx1_aln_ilm_{experiment}_{barcode}.bam.sort", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon1/medaka_{assembler}_pilon1.complete.flag", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.fasta", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/ilmn_fastqc_before/fastqc_complete.flag", experiment=EXPERIMENT, barcode=BARCODES),
        expand("results/{experiment}/{barcode}/dnadiff/medaka_{assembler}_pilon3.oriented.report", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk.bac120.summary.tsv", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_plasclass/plasclass.probs.out", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk_classification.txt", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/confirm_or_get_reference.flag", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/gtdbtk_sinfo_mod.json", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk_sinfo/update_sinfo_mod.flag", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/gtdbtk_sinfo_mod.json", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/medaka_{assembler}_dwnlds/download_done.flag", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/{assembler}_busco_ref/ref_busco_completion.flag", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_circos/{experiment}_{barcode}_circos.svg", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/assemblytics/{barcode}_medaka_{assembler}_pilon3.oriented.coords.csv", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/fastqc_nanopore/{experiment}_{barcode}_all_sqfilt_fastqc.html", experiment=EXPERIMENT, barcode=BARCODES),
        expand("results/{experiment}/{barcode}/nanoplot/{experiment}_{barcode}_all_sqfilt_NanoPlot-report.html", experiment=EXPERIMENT, barcode=BARCODES),
        expand("results/{experiment}/{barcode}/assembly_stats/medaka_{assembler}_pilon3.oriented.stats.txt", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),
        # expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_plasflow/medaka_{assembler}_pilon3_plasflow.tsv", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS),

rule preprocessing:
    input:
        "data/{experiment}/{barcode}"
    output:
        temp("results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq")
    log:
        "logs/{experiment}/{barcode}/preprocessing.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/preprocessing.py {input} {wildcards.experiment} {output} 2>&1 > {log}"


rule opt_read_filter:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    params:
        genSize = GENOME_SIZE,
        cov= COVERAGE
    output:
        temp("results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq")
    log:
        "logs/{experiment}/{barcode}/filter_optimization.log"
    conda:
        "envs/python3.yaml"
    shell:
        "python scripts/select_filter_params_by_cov_sm.py -f {input} -s {params.genSize} -c {params.cov} --dump fastq -decall 2>&1 > {log}"
        # "python scripts/select_filter_params_by_cov_sm.py -l 700 -q 9 -f {input} -s {params.genSize} -c 200 False --dump fastq 2>&1 > {log}" # for plasmid assembly
        # "python scripts/select_filter_params_by_cov_sm.py -l 700 -q 9 -f {input} -s {params.genSize} -opt False --dump fastq 2>&1 > {log}" # for plasmid assembly


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
        contigs="results/{experiment}/{barcode}/canu/{experiment}_{barcode}_canu.contigs.fasta",
        flag = "results/{experiment}/{barcode}/canu/canu.flag"
    params:
        genomeSize=GENOME_SIZE,
        useGrid="False"
    threads:
        16
    log:
        "logs/{experiment}/{barcode}/canu.log"
    conda:
        "envs/canu.yaml"
    shell:
        "canu -nanopore {input} -p {wildcards.experiment}_{wildcards.barcode}_canu -d {output.outdir} genomeSize={params.genomeSize} maxThreads={threads} useGrid={params.useGrid} 2>&1 > {log} && touch {output.flag}"


rule assemble_flye:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
    output:
        outdir=directory("results/{experiment}/{barcode}/flye"),
        outfile="results/{experiment}/{barcode}/flye/assembly.fasta",
        flag = "results/{experiment}/{barcode}/flye/flye.flag",
    params:
        genomeSize=230000,
        useGrid="False"
    threads:
        16
    log:
        "logs/{experiment}/{barcode}/canu.log"
    conda:
        "envs/flye.yaml"
    shell:
        "flye --nano-hq {input} -o {output.outdir} -g {params.genomeSize} -t {threads} 2>&1 > {log} && touch {output.flag}"


rule polish_medaka:
    input:
        reads="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta",
        asm_flag = "results/{experiment}/{barcode}/{assembler}/{assembler}.flag",
        draft_asm = get_draft_asm
    output:
        consensus="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta"
    params:
        outdir="results/{experiment}/{barcode}/medaka_{assembler}"
    log:
        "logs/{experiment}/{barcode}/medaka_{assembler}.log"
    threads:
        16
    conda:
        "envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.reads} -d {input.draft_asm} -o {params.outdir} -t {threads} 2>&1 > {log}"


rule dnadiff:
    input:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk.bac120.summary.tsv",
        asm="results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta",
        sinfo_flag="results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk_sinfo/update_sinfo_mod.flag",
        reference=get_reference_file
    output:
        "results/{experiment}/{barcode}/dnadiff/medaka_{assembler}_pilon3.oriented.report",
        "results/{experiment}/{barcode}/dnadiff/medaka_{assembler}_pilon3.oriented.delta"
    params:
        prefix="results/{experiment}/{barcode}/dnadiff/medaka_{assembler}_pilon3.oriented"
    conda:
        "envs/dnadiff.yaml"
    log:
        "logs/{experiment}/{barcode}/dnadiff/medaka_{assembler}_pilon3.oriented.log"
    shell:
        "dnadiff  -p {params.prefix} {input.reference} {input.asm} 2>&1 > {log}"


# rule prokka_expected_ref:
#     input:
#         # asm="results/{experiment}/{barcode}/medaka_{assembler}/consensus.fasta",
#         asm="results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta",
#         prot_file=get_ref_proteins
#     output:
#         protgbk_dir=directory("results/{experiment}/{barcode}/prokka_medaka_{assembler}_protgbk_exp"),
#         protg_tsv="results/{experiment}/{barcode}/prokka_medaka_{assembler}_protgbk_exp/{experiment}_{barcode}_{assembler}_prokka_protgbk_exp.tsv"
#     params:
#         genus=get_prokka_genus,
#         # species=get_prokka_species
#     log:
#         "logs/{experiment}/{barcode}/prokka_medaka_{assembler}_exp.log"
#     conda:
#         "envs/prokka.yaml"
#     shell:
#         # "prokka --kingdom Bacteria --genus {params.genus} --species {params.species} {input.asm} --usegenus --outdir {output.genusdb_dir} --prefix {wildcards.experiment}_{wildcards.barcode}_{wildcards.assembler}_prokka_genusdb --force && " 
#         "prokka --proteins {input.prot_file} {input.asm} --outdir {output.protgbk_dir} --prefix {wildcards.experiment}_{wildcards.barcode}_{wildcards.assembler}_prokka_protgbk_exp --force 2>&1 > {log}"


rule prokka_gtdbtk_ref:
    input:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk_classification.txt",
        "results/{experiment}/medaka_{assembler}_dwnlds/download_done.flag",
        asm="results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta",
        prot_file=get_ref_proteins
    output:
        protgbk_dir=directory("results/{experiment}/{barcode}/prokka_medaka_{assembler}_protgbk_class"),
        protg_tsv="results/{experiment}/{barcode}/prokka_medaka_{assembler}_protgbk_class/{experiment}_{barcode}_{assembler}_prokka_protgbk_class.tsv"
    log:
        "logs/{experiment}/{barcode}/prokka_medaka_{assembler}_class.log"
    conda:
        "envs/prokka.yaml"
    shell:
        # "prokka --kingdom Bacteria --genus {params.genus} --species {params.species} {input.asm} --usegenus --outdir {output.genusdb_dir} --prefix {wildcards.experiment}_{wildcards.barcode}_{wildcards.assembler}_prokka_genusdb --force && " 
        "prokka --proteins {input.prot_file} {input.asm} --outdir {output.protgbk_dir} --prefix {wildcards.experiment}_{wildcards.barcode}_{wildcards.assembler}_prokka_protgbk_class --force 2>&1 > {log}"


rule busco_assembly:
    input:
        asm="results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta"
    output:
        outfile="results/{experiment}/{barcode}/busco_medaka_{assembler}_pilon3/run_bacteria_odb10/short_summary.specific.{barcode}.txt"
    params:
        outdir="results/{experiment}/{barcode}/busco_medaka_{assembler}_pilon3",
        genus=BUSCO_GENUS,
        tmp_out="results/{experiment}/{barcode}/busco_medaka_{assembler}_pilon3/run_bacteria_odb10/short_summary.txt",
        summary_dir=get_busco_graph_outdir
    log:
        "logs/{experiment}/{barcode}/busco_medaka_{assembler}_pilon3.log"
    conda:
        "envs/busco.yaml"
    shell:
        "busco -i {input.asm} -l {params.genus} -o {params.outdir} -m genome -f && " # --auto-lineage-prok
        "mv {params.tmp_out} {output.outfile} 2>&1 > {log} && "
        "cp {output.outfile} {params.summary_dir}" 


rule busco_reference:
    input:
        expand("results/{experiment}/medaka_{assembler}_dwnlds/download_done.flag", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        references=get_references_path,
    params:
        genus=BUSCO_GENUS,
        summary_dir=get_busco_graph_outdir,
        outdir="results/{experiment}/{assembler}_busco_ref"
    output:
        "results/{experiment}/{assembler}_busco_ref/ref_busco_completion.flag"
    conda:
        "envs/busco.yaml"
    log:
        "logs/{experiment}/{assembler}_busco_ref/{assembler}_busco_reference.log"
    shell:
        "for ref_path in {input.references};"
        "do"
        " ref_name=$(basename $ref_path);"
        " IFS='.' read -ra SPLIT <<< $ref_name;"
        " ref_name=${{SPLIT[0]}};"
        " busco -i $ref_path -l {params.genus} -o {params.outdir}_$ref_name -m genome -f;"
        " mv {params.outdir}_$ref_name/run_bacteria_odb10/short_summary.txt {params.outdir}_$ref_name/run_bacteria_odb10/short_summary.specific.$ref_name.txt;"
        " cp {params.outdir}_$ref_name/run_bacteria_odb10/short_summary.specific.$ref_name.txt {params.summary_dir}; "
        "done; "
        "touch results/{wildcards.experiment}/{wildcards.assembler}_busco_ref/ref_busco_completion.flag 2>&1 > {log}"


rule busco_generate_plot:
    input:
        "results/{experiment}/{assembler}_busco_ref/ref_busco_completion.flag",
        expand("results/{experiment}/{barcode}/busco_medaka_{assembler}_pilon3/run_{genus}_odb10/short_summary.specific.{barcode}.txt", experiment=EXPERIMENT, assembler=ASSEMBLERS, genus=BUSCO_GENUS, barcode=BARCODES),
    output:
        "results/{experiment}/busco_graph/{assembler}/busco_figure.png",
    params:
        outdir = get_busco_graph_outdir
    conda:
        "envs/busco.yaml"
    shell:
        "python scripts/rename_files_for_busco_plot.py {params.outdir} && "
        "python scripts/busco_generate_plot.py -wd {params.outdir}" # 2>&1 > {log}"


rule fastqc_illumina_before:
    input:
        "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta", # this is not required as input for fastqc, but supplies the necessary wildcards and ensures 
        ilmn_reads=get_ilmn_reads
    output:
        "results/{experiment}/{barcode}/ilmn_fastqc_before/fastqc_complete.flag"
    params:
        outdir="results/{experiment}/{barcode}/ilmn_fastqc_before"
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/{experiment}/{barcode}/fastqc_ilmn_before_clipping.log"
    shell:
        "mkdir -p {params.outdir} && "
        "fastqc {input.ilmn_reads} -o {params.outdir} 2>&1 > {log} && "
        "touch {params.outdir}/fastqc_complete.flag"


rule fastqc_nanopore:
    input:
        ont_filt="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq",
        ont_raw="results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    output:
        "results/{experiment}/{barcode}/fastqc_nanopore/{experiment}_{barcode}_all_sqfilt_fastqc.html",
        "results/{experiment}/{barcode}/fastqc_nanopore/{experiment}_{barcode}_all_fastqc.html"
    params:
        outdir="results/{experiment}/{barcode}/fastqc_nanopore"
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/{experiment}/{barcode}/fastqc_nanopore.log"
    shell:
        "mkdir -p {params.outdir} && "
        "fastqc {input.ont_filt} -o {params.outdir} && "
        "fastqc {input.ont_raw} -o {params.outdir} 2>&1 > {log}"


rule nanoplot:
    input:
        ont_filt="results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fastq",
        ont_raw="results/{experiment}/{barcode}/{experiment}_{barcode}_all.fastq"
    output:
        "results/{experiment}/{barcode}/nanoplot/{experiment}_{barcode}_all_NanoPlot-report.html",
        "results/{experiment}/{barcode}/nanoplot/{experiment}_{barcode}_all_sqfilt_NanoPlot-report.html"
    params:
        outdir="results/{experiment}/{barcode}/nanoplot",
        pf_filt="{experiment}_{barcode}_all_sqfilt_",
        pf_raw="{experiment}_{barcode}_all_"
    conda:
        "envs/nanoplot.yaml"
    threads:
        4
    log:
        "logs/{experiment}/{barcode}/nanoplot.log"
    shell:
        "mkdir -p {params.outdir} && "
        "NanoPlot -t {threads} --fastq {input.ont_filt} -p {params.pf_filt} -o {params.outdir} && "
        "NanoPlot -t {threads} --fastq {input.ont_raw} -p {params.pf_raw} -o {params.outdir} 2>&1 > {log}"


checkpoint clip_adapters_sm:
    input:
        "results/{experiment}/{barcode}/ilmn_fastqc_before/fastqc_complete.flag", # this is not required as input for fastqc, but supplies the necessary wildcards and ensures 
        ilmn_reads=get_ilmn_reads
    output:
        directory("results/{experiment}/{barcode}/ilmn_clip_sm")
    params:
        adapters="resources/adapters_all_trimmomatic_condensed.fasta"
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/{experiment}/{barcode}/ilmn_clip_sm.log"
    shell:
        "mkdir -m 777 {output} -p && "
        "python3 scripts/identify_and_clip_ilmn_adapters_sm.py -cs 10 -a {params.adapters} -r {input.ilmn_reads} -o {output} 2>&1 > {log} && "
        "cd {output} && "
        "ls | grep _clipped.fq | xargs fastqc 2>&1 > fastqc_clip_adapters.log" #../../../../{log}" -o {output}


rule circlator_fixstart:
    input:
        "results/{experiment}/{barcode}/{assembly}.fasta"
    output:
        "results/{experiment}/{barcode}/circl_fixstart/{assembly}.oriented.fasta"#,
    params:
        out_prefix = "results/{experiment}/{barcode}/circl_fixstart/{assembly}.oriented"
    conda:
        "envs/circlator.yaml"
    log:
        "logs/{experiment}/{barcode}/{assembly}/circlator_fixstart.log"
    shell:
        "circlator fixstart {input} {params.out_prefix} 2>&1 > {log}"


rule assembly_stats:
    input:
        "results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta"
    output:
        "results/{experiment}/{barcode}/assembly_stats/medaka_{assembler}_pilon3.oriented.stats.txt"
    conda:
        "envs/python3.yaml"
    shell:
        "python3 scripts/assembly_statistics_4_bacterial.py {input} > {output}"


# rule circlator_all:
#     input:
#         asm = "results/{experiment}/{barcode}/{assembly}.fasta",
#         ont_reads = "results/{experiment}/{barcode}/{experiment}_{barcode}_all_sqfilt.fasta"
#     output:
#         "results/{experiment}/{barcode}/circl_all/{assembly}.oriented.fasta"#,
#         # directory("results/{experiment}/{barcode}/circlator")
#     params:
#         # out_prefix = "results/{experiment}/{barcode}/circl_all/{assembly}.oriented"
#         out_dir = directory("results/{experiment}/{barcode}/circl_all")
#     conda:
#         "envs/circlator.yaml"
#     log:
#         "logs/{experiment}/{barcode}/{assembly}/circlator_all.log"
#     shell:
#         "circlator all {input.asm} {input.ont_reads} {params.out_dir} 2>&1 > {log}"


rule classify_gtdbtk:
    input:
        asm = "results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta"
    output:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk.bac120.summary.tsv"
    params:
        genome_dir = "results/{experiment}/{barcode}/medaka_{assembler}_pilon3",
        out_dir= "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk"
    threads:
        8
    conda:
        "envs/gtdbtk.yaml"
    log:
        "logs/{experiment}/{barcode}/medaka_{assembler}_pilon3/classify_gtdbtk.log"
    shell:
        # "gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir> 2>&1 > {log}"
        "mamba env config vars set GTDBTK_DATA_PATH=/home/simon/mambaforge/envs/gtdbtk/share/gtdbtk-2.3.2/db && "
        # "mamba env config vars set GTDBTK_DATA_PATH=/homes/simon/.conda/envs/gtdbtk/share/gtdbtk-2.3.2/db && "
        "gtdbtk classify_wf --genome_dir {params.genome_dir} --out_dir {params.out_dir} --extension fasta --cpu {threads} --mash_db resources 2>&1 > {log}"


rule set_gtdbtk_classification:
    input:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk.bac120.summary.tsv"
    output:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk_classification.txt",
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_update-sinfo/Sample_Info_{barcode}.json"
    run:
        calls = []
        with open(input[0], "r") as gtdbtk_out, open(output[0], "w") as out:
            data = gtdbtk_out.read().split("\n")
            keys = data[0].strip().split("\t")
            data = [l for l in data if len(l.strip())>0][1:]
            # print(f"INFO: {input} contains {len(data)} call(s)")
            for line in data:
                data = line.strip().split("\t")
                ddct = dict(zip(keys, data))
                calls.append(ddct)
                class_tree = ddct["classification"]
                class_dct = dict(
                    [call.split("__") for call in class_tree.split(";")]
                    )
                # print("class_dct", class_dct)
                class_dct = {GTDB_TRANSL[k]:v for k,v in class_dct.items()}
                # print("class_dct", class_dct)
                class_cats = [
                    "gtdb_species",
                    "gtdb_genus",
                    "gtdb_family",
                    "gtdb_order",
                    "gtdb_class",
                    "gtdb_phylum",
                    "gtdb_domain",
                    ]
                spec_call = False
                for cat in class_cats:
                    try:
                        gtdbtk_call = class_dct[cat]
                        # print(f"INFO: Using {cat} as gtbdtk call category.")
                        # print(f"GTDBTK call: {class_dct[cat]}")
                        if cat=="gtdb_species":
                            spec_call = True
                        break
                    except:
                        pass
                        # print(f"WARNING: gtdbtk did not yield a {cat} call.")
                out.write(gtdbtk_call)
                for cat in class_dct:
                    SAMPLE_INFO[wildcards.barcode][cat] = class_dct[cat]
                SAMPLE_INFO[wildcards.barcode]["gtdb_ref"] = ddct["fastani_reference"]
                # print(f"SAMPLE_INFO: {SAMPLE_INFO}")
                # print(
                #     f"Suggested reference: {SAMPLE_INFO[wildcards.barcode]['gtdb_ref']}"
                #     )
                with open(output[1], "w") as sinfo:
                    augm_sinfo = json.dumps(
                        {wildcards.barcode: SAMPLE_INFO[wildcards.barcode]},
                        indent=4
                        )
                    sinfo.write(augm_sinfo)


checkpoint update_sample_info_gtdbtk_class:
    input:
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_update-sinfo/Sample_Info_{barcode}.json", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS)
    output:
        expand("results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/gtdbtk_sinfo_mod.json", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        expand("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk_sinfo/update_sinfo_mod.flag", experiment=EXPERIMENT, barcode=BARCODES, assembler=ASSEMBLERS)
    run:
        # print("Input of update_sample_info_gtdbtk_class", input)
        # print("output of update_sample_info_gtdbtk_class", output)

        inf_by_asm = defaultdict(list)
        for asm in ASSEMBLERS:
            for in_file in input:
                if asm in in_file:
                    inf_by_asm[asm].append(in_file)

        for asm in inf_by_asm:
            for in_file in inf_by_asm[asm]:
                with open(in_file, "r") as inf:
                    bc_sinfo = json.loads(inf.read())
                    for key in bc_sinfo:
                        SAMPLE_INFO[key] = bc_sinfo[key]
                if_spl = in_file.split("/")
                ofs = if_spl.index("results")
                exp, bc, asm =  if_spl[ofs+1], if_spl[ofs+2], if_spl[ofs+3].split("_")[1]
                flag = f"results/{exp}/{bc}/medaka_{asm}_pilon3_gtdbtk_sinfo/update_sinfo_mod.flag"
                with open(flag, "w") as sflag:
                    sflag.write(flag)
            # print(f"SAMPLE_INFO for {asm} from update_sample_info_gtdbtk_class", SAMPLE_INFO)

            if not os.path.exists(f"results/{exp}/medaka_{asm}_pilon3_gtdbtk_sinfo"):
                os.mkdir(f"results/{exp}/medaka_{asm}_pilon3_gtdbtk_sinfo")
            with open(os.path.join("results", exp, f"medaka_{asm}_pilon3_gtdbtk_sinfo", "gtdbtk_sinfo_mod.json"), "w") as sinfo:
                sinfo_mod = json.dumps(SAMPLE_INFO, indent=4)
                sinfo.write(sinfo_mod)


checkpoint confirm_or_get_reference:
    input:
        "results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/gtdbtk_sinfo_mod.json"
    output:
        "results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/confirm_or_get_reference.flag",
        dwnl_dir=directory("results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_dwnld")
    run:
        with open(
            os.path.join(
                "results",
                wildcards.experiment,
                f"medaka_{wildcards.assembler}_pilon3_gtdbtk_sinfo",
                "gtdbtk_sinfo_mod.json"
                ), 
            "r"
            ) as sinfo:
            sinfo_mod = json.load(sinfo)
            # print("sinfo_mod", sinfo_mod)
        dwnld = set()
        for sample in sinfo_mod:
            messages = []
            best_ref = sinfo_mod[sample]["gtdb_ref"]
            ref_dir = os.path.join("resources", EXPERIMENT)
            avail_ref = glob.glob(os.path.join(ref_dir, f"*{best_ref}*"))
            req_ext = ["fna", "fasta", "fa", "faa", "gbk", "gb"]
            if len(avail_ref) > 0:
                m1 = f"Signature of best ref for {sample} assembled with {wildcards.assembler} already present in resources"
                messages.append(m1)
                # print(m1) 
            corr_ext = [f.split(".")[-1] for f in avail_ref if f.split(".")[-1] in req_ext]
            ref_type = {"fna": "asm", "fasta": "asm", "fa": "asm", "faa": "prot", "gbk": "prot", "gb": "prot"}
            ref_set = {ref_type[ext] for ext in corr_ext}
            if len(ref_set) >= 2:
                m2 = "It appears all required reference files are present."
            else:
                m2 = f"Some reference file/s is/are missing. Adding {best_ref} to download set."
                dwnld.add(best_ref)
            messages.append(m2)
            # print(m2)
            out_flag = open(output[0], "w")
            out_flag.close()
            with open(f"results/{EXPERIMENT}/medaka_{wildcards.assembler}_pilon3_gtdbtk_sinfo/confirm_or_get_reference.flag", "a") as flag:
                for m in messages:
                    flag.write(f"{m}\n")
        # print("dwnld_set", dwnld)
        os.makedirs(f"results/{EXPERIMENT}/medaka_{wildcards.assembler}_pilon3_gtdbtk_dwnld", exist_ok=True)
        for ref in dwnld:
            # print("Writing download file for: ", ref)
            with open(f"results/{EXPERIMENT}/medaka_{wildcards.assembler}_pilon3_gtdbtk_dwnld/{ref}.dwnld", "w") as lst:
                lst.write(ref)


checkpoint download_references:
    input:
        "results/{experiment}/medaka_{assembler}_pilon3_gtdbtk_sinfo/confirm_or_get_reference.flag",
        accs=get_genome_dwnld_files
    output:
        "results/{experiment}/medaka_{assembler}_dwnlds/download_done.flag",
    conda:
        "envs/datasets.yaml"
    shell:
        "for acc in {input.accs};"
        "do"
        " ref_acc=$(grep _ $acc);"
        " mkdir -p resources/{wildcards.experiment}/$ref_acc;"
        # do not include cds download, because the cds file end in fna, confounding detection of assemblies
        " datasets download genome accession $ref_acc --include genome,protein,gff3,rna,gbff,gtf,seq-report --filename resources/{wildcards.experiment}/$ref_acc.zip;"
        " unzip -n resources/{wildcards.experiment}/$ref_acc.zip -d resources/{wildcards.experiment}/$ref_acc;"
        # The inner one of the nested for loops was adopted with changes from
        # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/rename-files/
        " for file in resources/{wildcards.experiment}/$ref_acc/ncbi_dataset/data/*/*;"
            "do"
            " directory_name=$(dirname $file);"
            " accession=$(basename $directory_name);"
            " mv ${{file}} ${{directory_name}}/${{accession}}_$(basename $file);"
            "done;"
        " mv -n resources/{wildcards.experiment}/$ref_acc/ncbi_dataset/data/$ref_acc/*  resources/{wildcards.experiment};"
        "done && "
        "touch results/{wildcards.experiment}/medaka_{wildcards.assembler}_dwnlds/download_done.flag"


# rule plasflow:
#     input:
#         asm = "results/{experiment}/{barcode}/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.fasta"
#     output:
#         # "results/{experiment}/{barcode}/circl_fixstart/{assembly}.oriented.fasta"#,
#         # directory("results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk")
#         "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_plasflow/medaka_{assembler}_pilon3_plasflow.tsv"
#     params:
#         # out_prefix = "results/{experiment}/{barcode}/circl_fixstart/{assembly}.oriented"
#         # genome_dir = "results/{experiment}/{barcode}/medaka_{assembler}_pilon3",
#         # out_dir= "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk"
#     conda:
#         "envs/plasflow.yaml"
#     log:
#         "logs/{experiment}/{barcode}/medaka_{assembler}_pilon3/plasflow.log"
#     shell:
#         # "curr_env=$(conda info --env | grep '*' | rev | cut -d' ' -f1 | rev) && "
#         # "python $curr_env/lib/python3.5/site-packages/PlasFlow/PlasFlow.py --input {input} --output {output} 2>&1 > {log}"
#         "python pkgs/PlasFlow/PlasFlow.py --input {input} --output {output} 2>&1 > {log}"


rule plasclass:
    input:
        asm = "results/{experiment}/{barcode}/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.fasta" # consider using .oriented
    output:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_plasclass/plasclass.probs.out"
    threads:
        8
    conda:
        "envs/plasclass.yaml"
    log:
        "logs/{experiment}/{barcode}/medaka_{assembler}_pilon3/plasclass.log"
    shell:
        # "curr_env=$(conda info --env | grep '*' | rev | cut -d' ' -f1 | rev) && "
        # "python $curr_env/lib/python3.9/site-packages/plasclass/plasclass.py -f {input} -o {output} -p {threads}"
        "python pkgs/PlasClass/classify_fasta.py -f {input} -o {output} -p {threads}"
        

rule mummer2circos_single:
    input:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk/gtdbtk.bac120.summary.tsv",
        sinfo_flag="results/{experiment}/{barcode}/medaka_{assembler}_pilon3_gtdbtk_sinfo/update_sinfo_mod.flag",
        dwnld_flag=expand("results/{experiment}/medaka_{assembler}_dwnlds/download_done.flag", experiment=EXPERIMENT, assembler=ASSEMBLERS),
        asm="results/{experiment}/{barcode}/circl_fixstart/medaka_{assembler}_pilon3/medaka_{assembler}_pilon3.oriented.fasta",
        ref=get_reference_file
    output:
        "results/{experiment}/{barcode}/medaka_{assembler}_pilon3_circos/{experiment}_{barcode}_circos.svg"
    params:
        prefix="{experiment}_{barcode}_circos",
        workdir="results/{experiment}/{barcode}/medaka_{assembler}_pilon3_circos"
    threads:
        4
    conda:
        # "envs/mummer2circos_pkg.yaml"
        "envs/mummer2circos.yaml"
    log:
        "logs/{experiment}/{barcode}/medaka_{assembler}_pilon3/mummer2circos.log"
    shell:
        "mkdir -p {params.workdir} && "
        "cd {params.workdir} && "
        "mummer2circos -f -l -r ../../../../{input.ref} -q ../../../../{input.asm} -o {params.prefix}"
        #"python3 pkgs/mummer2circos/mummer2circos/mummer2circos.py -l -r {input.ref} -q {input.asm} -o {params.prefix}"


rule assemblytics:
    input:
        "results/{experiment}/{barcode}/dnadiff/medaka_{assembler}_pilon3.oriented.delta"
    output:
        # "results/{experiment}/{barcode}/assemblytics/{barcode}_medaka_{assembler}_pilon3.oriented.Assemblytics.Dotplot_filtered.png"
        "results/{experiment}/{barcode}/assemblytics/{barcode}_medaka_{assembler}_pilon3.oriented.coords.csv"
    params:
        output_prefix="results/{experiment}/{barcode}/assemblytics/{barcode}_medaka_{assembler}_pilon3.oriented",
        unique_anchor_length=10000,
        min_variant_size=10000,
        max_variant_size=1,
    conda:
        "envs/assemblytics.yaml"
    shell:
        "Assemblytics {input} {params.output_prefix} {params.unique_anchor_length} {params.min_variant_size} {params.max_variant_size}"

# add a third round of pilon polishing; checked
# classification; checked
# Plasmids; checked
# transfer BUSCO plot generation downstream of classifciation; checked
# Synteny; checked - covered by assemblytics
# Detect and analyze variants: assemblytics; checked
# Circular plot with GC skew etc.; checked
# Quality assessment (QUAST or similar); probably mostly covered by assemblytics
# find a way to remove clipped illumina reads (temp() not working due to unknown name)
# write a rule to generate assembly graphs from flye (canu also?) assembly_info.gv files using graphviz with following commands:
    # dot {results}/{barcode}/flye/assembly_graph.gv > {results}/{barcode}/flye/{barcode}_asm_graph.dot &&
    # dot -Tpng {results}/{barcode}/flye/{barcode}_asm_graph.dot -o {results}/{barcode}/flye/{barcode}_asm_graph.png
# use {results}/{barcode}/flye/assembly_graph.gv as input file, which would necessitate to add it to flye output files as graph=...
# migth also be tried with gfa file
# alternative for graphviz could be bandage

