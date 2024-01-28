#! /usr/local/bin/python3

import sys
from collections import defaultdict, Counter
import statistics as stat
from matplotlib import pyplot as plt


# # original from intloc.il_aux.il_integrator func: parse_acceptor_seq(acceptor)
# def parse_assembly(fasta):
#     """Scan assembly into which the integrations will be inserted"""
#     contigs = dict()
#     curr_cont = ""
#     n_stretch = False
#     curr_ns_start = None
#     n_regions = defaultdict(list)
#     with open(fasta, "r") as acc:
#         for line in acc:
#             if line.startswith(">"):
#                 if n_stretch:
#                     n_stretch = False
#                     cns_end = contigs[curr_cont]
#                     n_regions[curr_cont].append((curr_ns_start, cns_end))
#                 curr_cont = line.strip()
#                 contigs[curr_cont] = 0
#             else:
#                 if n_stretch:
#                     bases = ["A","T","G","C","a","t","g","c"]
#                     chk_end = [base for base in line if base in bases]
#                     if any(chk_end):
#                         n_stretch = False
#                         cns_end = min([line.index(base) for base in chk_end])
#                         cns_end = cns_end + contigs[curr_cont]
#                         n_regions[curr_cont].append((curr_ns_start, cns_end))
#                 else:
#                     if "N" in line:
#                         n_stretch = True
#                         curr_ns_start = contigs[curr_cont] + line.index("N")
#                 contigs[curr_cont] += len(line.strip())
#         if n_stretch:
#                     n_stretch = False
#                     cns_end = contigs[curr_cont]
#                     n_regions[curr_cont].append((curr_ns_start, cns_end))

#     stats = [contigs[tig] for tig in contigs]
#     total_bp = sum(stats)
#     gaps = [reg for reg_lst in n_regions.values() for reg in reg_lst]
#     gap_lens = [reg[1]-reg[0] for reg in gaps]
#     total_ns = sum([reg[1]-reg[0] for reg in gaps])

#     # print(gap_lens)
#     # print(contigs)

#     print("Acceptor Seq. statistics:")
#     print("No. of contigs", len(stats))
#     print("Average contig length:", stat.mean(stats))
#     print("longest contig", max(stats))
#     print("Total base count", total_bp)
#     print("Number of gaps (N-stretches) in genome", len(gaps))
#     print("Total N count in genome", total_ns)
#     if len(gap_lens) > 0:
#         print("Average gap length", sum(gap_lens)/len(gap_lens))

#     return contigs, total_bp, n_regions

def analyze_tig(seq):
    len_bp = len(seq)
    base_comp = Counter(seq)
    idx = -1
    n_stretch = False
    curr_ns = []
    n_stretch_pos = []
    for char in seq:
        idx += 1
        if char == "N":
            if n_stretch:
                pass
            else:
                n_stretch = True
                curr_ns.append(idx)
        else:
            if n_stretch:
                n_stretch = False
                curr_ns.append(idx)
                n_stretch_pos.append(curr_ns)
                curr_ns = []
            else:
                pass
    if n_stretch:
        n_stretch = False
        curr_ns.append(idx + 1)
        n_stretch_pos.append(curr_ns)
        curr_ns = []
    return len_bp, base_comp, n_stretch_pos


def plot_readlen(read_lens: list, fasta_in: str):
    len_ct = Counter(read_lens)
    len_ct = {k:v for k,v in len_ct.items() if v>1}
    print(f"len_count: {len_ct}")
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.bar(list(len_ct.keys()), list(len_ct.values()))
    if len(len_ct) > 0:
        min_len, max_len = min(list(len_ct.keys())), max(list(len_ct.keys()))+1
        tick_range = [tick for tick in range(min_len, max_len)]
    # plt.xticks(tick_range, tick_range, rotation=90)
    # ax.set_xticklabels(tick_range, rotation=90)
    plt.tight_layout()
    plt.savefig(f"{fasta_in}_readlen_plot.png")


def parse_assembly(fasta):
    """Scan assembly into which the integrations will be inserted"""
    # contigs = {}
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
    all_bases = Counter()
    all_N_runs = []
    
    for tig in cont_info:
        all_lens.append(cont_info[tig][0])
        all_bases.update(cont_info[tig][1])
        all_N_runs.extend(cont_info[tig][2])
        # print(tig)
        # print(
        #     f"len_bp: {cont_info[tig][0]},",
        #     f"base_comp: {cont_info[tig][1]},",
        #     f"n_stretch_pos: {cont_info[tig][2]}"
        #     )

    print(all_lens)
    plot_readlen(all_lens, fasta)
    total_bp = sum(all_lens)
    gap_lens = [reg[1]-reg[0] for reg in all_N_runs]

    print("Assembly statistics:")
    print("No. of contigs", len(all_lens))
    print("Average contig length:", stat.mean(all_lens))
    print("longest contig", max(all_lens))
    print("Total base count", total_bp)
    print("Number of gaps (N-stretches) in genome", len(all_N_runs))
    print("Total N count in genome", all_bases["N"])
    if len(gap_lens) > 0:
        print("Average gap length", int(sum(gap_lens)/len(gap_lens)))

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