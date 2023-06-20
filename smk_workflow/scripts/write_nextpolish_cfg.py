

mtj = "{multithread_jobs}"

cfg_templ = f"[General]\n"\
        f"job_type = local\n"\
        f"job_prefix = nextPolish\n"\
        f"task = best\n"\
        f"rewrite = yes\n"\
        f"rerun = 3\n"\
        f"parallel_jobs = 6\n"\
        f"multithread_jobs = 5\n"\
        f"genome = ./raw.genome.fasta #genome file\n"\
        f"genome_size = auto\n"\
        f"workdir = ./01_rundir\n"\
        f"polish_options = -p {mtj}\n"\
        f"\n"\
        f"[lgs_option]\n"\
        f"lgs_fofn = ./lgs.fofn\n"\
        f"lgs_options = -min_read_len 1k -max_depth 100\n"\
        f"lgs_minimap2_options = -x map-ont"

print(cfg_templ)
