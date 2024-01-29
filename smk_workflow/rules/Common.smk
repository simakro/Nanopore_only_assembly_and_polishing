
def get_ref_path(wildcards):
    ref_name = SAMPLE_INFO[wildcards.barcode]["ref"]
    return os.path.join("resources", wildcards.experiment, ref_name)


def get_ref_proteins(wildcards): # adapt to also accept .gb (which is same format as gbk) extension and .faa which is a different pure protein format usable for prokka --proteins
    reference = SAMPLE_INFO[wildcards.barcode]["ref"]
    comp_ext = [".gbk", ".gb", ".faa"]
    prot_files = [f'{(".").join(reference.split(".")[:-1])}{ext}' for ext in comp_ext]
    res_dir = os.path.join("resources", wildcards.experiment)
    res_files = [entry.name for entry in os.scandir(res_dir)]
    avail_prot = [f for f in prot_files if f in res_files][0]
    # ref_path = os.path.join("resources", wildcards.experiment, SAMPLE_INFO[wildcards.barcode]["ref"])
    prot_path = os.path.join(res_dir, avail_prot)
    return prot_path


def get_draft_asm(wildcards):
    if wildcards.assembler == "canu":
        return f"results/{wildcards.experiment}/{wildcards.barcode}/canu/{wildcards.experiment}_{wildcards.barcode}_canu.contigs.fasta"
    elif wildcards.assembler == "flye":
        return f"results/{wildcards.experiment}/{wildcards.barcode}/flye/assembly.fasta"
    else:
        print("WARNING: No assembler in wildcards. Defaulting to canu.")
        return "canu"


def get_prokka_genus(wildcards):
    return SAMPLE_INFO[wildcards.barcode]["genus_fallback"]


def get_prokka_species(wildcards):
    return SAMPLE_INFO[wildcards.barcode]["species_fallback"]


def get_len_filt_param(wildcards):
    param_json = f"results/{wildcards.experiment}/{wildcards.barcode}/filt_params.json"
    with open(param_json, "r") as paraj:
        data = paraj.read()
        param_dct = json.loads(data)
        print(param_dct.items())
    return param_dct["len"]


def get_qual_filt_param(wildcards):
    param_json = f"results/{wildcards.experiment}/{wildcards.barcode}/filt_params.json"
    with open(param_json, "r") as paraj:
        data = paraj.read()
        param_dct = json.loads(data)
        print(param_dct.items())
    return param_dct["qual"]


def get_busco_graph_outdir(wildcards):
    # return os.path.split(wildcards.input)[0]
    # outdir = f"results/{wildcards.experiment}/busco_graph_{wildcards.assembler}"
    outdir = f"results/{wildcards.experiment}/busco_graph"
    if not os.path.exists(f"results/{wildcards.experiment}"):
       os.mkdir(f"results/{wildcards.experiment}")
    if not os.path.exists(outdir):
       os.mkdir(outdir)
    print(outdir)
    return outdir


def get_ilmn_reads(wildcards):
    files = SAMPLE_INFO[wildcards.barcode]["illumina"]
    path = os.path.join("data", wildcards.experiment, "short-reads")
    paths = [os.path.join(path, f) for f in files]
    # both_space_sep = f"{paths[0]} {paths[1]}" # surprisingly returning the paths-list appears to work as input for all cli commands; possibly snakemake already handles converting them in space seperated strings?
    return paths


def get_clipped_ilmn_reads(wildcards):
    ilmn_clip_dir = checkpoints.clip_adapters_sm.get(**wildcards).output[0]
    # ilmn_reads = get_ilmn_reads(wildcards)
    return glob.glob(os.path.join(ilmn_clip_dir, "*_clipped.fq"))