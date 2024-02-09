#! /usr/local/bin/python3

import sys
from collections import defaultdict, Counter
import statistics as stat


def analyze_tig(seq):
    len_bp = len(seq)
    base_comp = Counter(seq.upper())
    gc_content = ((base_comp["G"]+base_comp["C"])/len_bp)*100
    return len_bp, gc_content


def parse_assembly(fasta):
    """Scan assembly into which the integrations will be inserted"""
    curr_name = ""
    curr_seq = ""
    cont_info = {}
    n_regions = defaultdict(list)
    with open(fasta, "r") as acc:
        for line in acc:
            if line.startswith(">"):
                if len(curr_name) > 0:
                    cont_info[curr_name] = analyze_tig(curr_seq)
                curr_name = line.strip()
                curr_seq = ""
            else:
                curr_seq += line.strip()
        cont_info[curr_name] = analyze_tig(curr_seq)

    all_lens = []
    # all_bases = Counter()
    
    print("\nContig stats:")
    print("name,length,GC%")
    for tig in cont_info:
        print(tig[1:], cont_info[tig][0], cont_info[tig][1])
        all_lens.append(cont_info[tig][0])
        # all_bases.update(cont_info[tig][1])
    print("\n")
    total_bp = sum(all_lens)

    print("Assembly statistics:")
    print("No. of contigs", len(all_lens))
    print("Average contig length:", stat.mean(all_lens))
    print("longest contig", max(all_lens))
    print("Total base count", total_bp)
    return cont_info, total_bp, n_regions


def calculate_Nx(x, cont_info, total_bp):
    cont_lens = [cont_info[tig][0] for tig in cont_info]
    threshold = total_bp * (x/100)
    tigs = sorted(cont_lens)
    tigs = tigs[::-1]
    curr_val = 0
    curr_tig = -1
    while curr_val <= threshold:
        curr_tig += 1
        curr_val += tigs[curr_tig]
    Nx = tigs[curr_tig]
    return Nx


if __name__ == "__main__":
    fasta = sys.argv[1]
    contigs, total_bp, n_regions = parse_assembly(fasta)
    x = 50
    Nx = calculate_Nx(x, contigs, total_bp)
    print(f"N{x} = {Nx}")