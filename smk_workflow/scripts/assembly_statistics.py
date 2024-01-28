#! /usr/local/bin/python3

import sys
import statistics as stat

def parse_assembly(fasta):
    contigs = {}
    curr_cont = ""
    with open(fasta, "r") as acc:
        for line in acc:
            if line.startswith(">"):
                curr_cont = line.strip()
                contigs[curr_cont] = 0
            else:
                contigs[curr_cont] += len(line.strip())
    stats = [contigs[tig] for tig in contigs]
    total_bp = sum(stats)
    #print(contigs)
    for tig in contigs:
        print(f"{tig}:\n{contigs[tig]}")
    print("Assembly statistics:")
    print("No. of contigs", len(stats))
    print("Average contig length:", stat.mean(stats))
    print("longest contig", max(stats))
    print("Total base count", total_bp)
    return contigs, total_bp

def calculate_Nx(x, contigs, total_bp):
    threshold = total_bp * (x/100)
    tigs = sorted(contigs.values())
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
    contigs, total_bp = parse_assembly(fasta)
    x = 50
    Nx = calculate_Nx(x, contigs, total_bp)
    print(f"N{x} = {Nx}")

