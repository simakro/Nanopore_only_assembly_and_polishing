
def get_references_path(wildcards):
    print("Getting references path")
    ref_path=f"resources/{wildcards.experiment}"
    references = []
    fasta_exts = ["fa", "fna", "fasta"]
    for ext in fasta_exts:
        ext_refs=glob.glob(os.path.join(ref_path, f"*.{ext}"))
        references.extend(ext_refs)
    print("references", references)
    return references


def get_reference_file(wildcards):
    print("Getting reference file")
    sinfo_mod = f"results/{wildcards.experiment}/medaka_{wildcards.assembler}_pilon2_gtdbtk_sinfo/gtdbtk_sinfo_mod.json"
    if os.path.exists(sinfo_mod):
        print("Updated sample_info in json format is available")
        try:
            print("Trying to retrieve reference based on gtdbtk classification")
            ref_acc = SAMPLE_INFO[wildcards.barcode]["gtdb_ref"]
            print("ref_acc", ref_acc)
            ref_path=glob.glob(os.path.join("resources", wildcards.experiment, f"{ref_acc}_*_genomic.fna"))[0]
            print(f"Reference file path: {ref_path}")
        except:
            print("reference file path could not be determined. Reverting to fallback.")
            ref_name = SAMPLE_INFO[wildcards.barcode]["ref_fallback"]
            ref_path = os.path.join("resources", wildcards.experiment, ref_name)
            
    else:
        print("Updated sample_info in json format does not yet exist")
        print("Using fallback for the time being")
        ref_name = SAMPLE_INFO[wildcards.barcode]["ref_fallback"]
        ref_path = os.path.join("resources", wildcards.experiment, ref_name)
    return ref_path


def get_fallback_ref(wildcards):
    reference = SAMPLE_INFO[wildcards.barcode]["ref"]
    comp_ext = [".gbk", ".gb", ".faa", ".gbff"]
    prot_files = [f'{(".").join(reference.split(".")[:-1])}{ext}' for ext in comp_ext]
    res_dir = os.path.join("resources", wildcards.experiment)
    res_files = [entry.name for entry in os.scandir(res_dir)]
    try:
        avail_prot = [f for f in prot_files if f in res_files][0]
        prot_path = os.path.join(res_dir, avail_prot)
    except:
        prot_path = f"results/{wildcards.experiment}/medaka_{wildcards.assembler}_dwnlds/download_done.flag"
    return prot_path


def get_ref_proteins(wildcards):  
    sinfo_mod = f"results/{wildcards.experiment}/medaka_{wildcards.assembler}_pilon2_gtdbtk_sinfo/gtdbtk_sinfo_mod.json"
    print("sinfo_mod", sinfo_mod)
    if os.path.exists(sinfo_mod):
        print("Updated sample_info in json format exists")
        with open(sinfo_mod, "r") as si_mod:
            Sample_Info_mod = json.loads(si_mod.read())
        gtdb_ref = Sample_Info_mod[wildcards.barcode]["gtdb_ref"]
        gbfl=glob.glob((os.path.join("resources", wildcards.experiment, f"{gtdb_ref}_genomic.gbff")))
        if len(gbfl)>0:
            print("GBFF file for reference found in resources")
            prot_file = gbfl[0]
        else:
            print("No GBFF file for reference found. Using fallback.")
            prot_file = get_fallback_ref(wildcards)
    else:
        print("Updated sample_info-json does not exist. Using fallback.")
        prot_file = get_fallback_ref(wildcards)
    return prot_file


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
    outdir = f"results/{wildcards.experiment}/busco_graph/{wildcards.assembler}"
    if not os.path.exists(f"results/{wildcards.experiment}"):
       os.mkdir(f"results/{wildcards.experiment}")
    if not os.path.exists(f"results/{wildcards.experiment}/busco_graph"):
       os.mkdir(f"results/{wildcards.experiment}/busco_graph")
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


def get_genome_dwnld_files(wildcards):
    dwnld_dir = checkpoints.confirm_or_get_reference.get(**wildcards).output.dwnl_dir
    print("dwnld_dir", dwnld_dir)
    # ilmn_reads = get_ilmn_reads(wildcards)
    dwnl_paths = glob.glob(os.path.join(dwnld_dir, "*.dwnld"))
    print("dwnl_paths", dwnl_paths)
    # dwnl_fls = [os.path.split(f)[-1] for f in dwnl_paths]
    # print("dwnl_fls", dwnl_fls)
    # accessions = [".".join(f.split(".")[:-1]) for f in dwnl_fls]
    # print("accessions", accessions)
    return dwnl_paths