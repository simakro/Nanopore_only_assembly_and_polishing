import os
import sys


def construct_filenames(kraken_outdir):
    class_file = os.path.join(kraken_outdir, "classifications_nonhost_reads.kraken")
    report = os.path.join(kraken_outdir, "classifications_nonhost_reads.report")
    return report, class_file


def read_kraken_file(k_file: str) -> dict:
    kraken_class = {"C": {}, "U": {}}
    header = ["tax-id", "seq_len"]
    with open(k_file, "r") as kf:
        for line in kf:
            l = line.strip()
            if len(l):
                ls = l.split("\t")
                class_status, seq_id = ls[0], ls[1]
                info = dict(zip(header, ls[2:4]))
                # info = dict(zip(header, [int(n) for n in ls[2:4]]))
                info["seq_len"] = int(info["seq_len"])
                kraken_class[class_status][seq_id] = info
    # print(kraken_class)
    return kraken_class


def read_kraken_report(report: str) -> dict:
    ranks = {
        "R": "Root",
        "D": "Domain",
        "K": "Kingdom",
        "P": "Phylum",
        "C": "Class",
        "O": "Order",
        "F": "Family",
        "G": "Genus",
        "S": "Species",
    }
    header = [
        "percent_reads",
        "cum_num_reads",
        "num_class_reads",
        "taxon_rank",
        "ncbi-tax-id",
        "scientific_name"
    ]
    with open(report, "r") as rep:
        all_lines = []
        for line in rep:
            line = line.strip()
            if line:
                ls = line.split("\t")
                ls = [col.strip() for col in ls]
                ld = dict(zip(header, ls))
                all_lines.append(ld)
    final_class_lines = [l for l in all_lines if int(l["num_class_reads"]) > 0]
    # print(final_class_lines)
    return final_class_lines


def associate_taxnames_with_readnames(
    report_info: dict, class_file_data: dict, sort_key: str = "seq_len"  # "contig_id"
):
    report = []
    classified_tigs = class_file_data["C"]
    # filter out human read using ncbi tax-id for homo sapiens
    non_human_reads = {k:v for k,v in classified_tigs.items() if int(v["tax-id"]) != 9606}
    # filter out any reads that are not classified at least at family level
    req_rank = ["S", "G", "F"]
    filt_report_info = [r for r in report_info if r["taxon_rank"][0] in req_rank]
    excl_for_rank = len(report_info) - len(filt_report_info)
    ct_human_reads = len(classified_tigs) - len(non_human_reads)
    unclassified = class_file_data["U"]
    for tig in non_human_reads:
        tig_id = tig
        tig_info = classified_tigs[tig]
        tig_info["contig_id"] = tig_id
        for class_line in filt_report_info:
            if class_line["ncbi-tax-id"] == tig_info["tax-id"]:
                tig_info["class_rank"] = class_line["taxon_rank"]
                tig_info["class_name"] = class_line["scientific_name"]
                report.append(tig_info)
    report = sorted(report, key=lambda x: x[sort_key])[::-1]
    # additional info messages
    unclass_msg = f"# {len(unclassified)} reads could not be classified"
    filt_hs_reads_msg = f"# {ct_human_reads} human reads not included in this summary"
    rank_excl_msg = f"# {excl_for_rank} reads not classified at least at family level were excluded"
    msgs = [unclass_msg, filt_hs_reads_msg, rank_excl_msg]
    report.append("")
    for m in msgs:
        if len(m):
            report.append(m)
    # print("report", report)
    return report


def write_custom_report(report_data: list, report_file: str):
    header = "contig\tseq_len\tclass_rank\tclassification"
    with open(report_file, "w") as rep:
        rep.write(header + "\n")
        for tig in report_data:
            if type(tig) == dict:
                l = f'{tig["contig_id"]}\t{tig["seq_len"]}\t{tig["class_rank"]}\t{tig["class_name"]}'
                rep.write(l + "\n")
            elif type(tig) == str:
                rep.write(tig + "\n")
            else:
                print("Encountered unexpected type in classification report data")


if __name__ == "__main__":
    kraken_outdir = sys.argv[1]
    post_proc_report = os.path.join(kraken_outdir, "custom_summary.tsv") 
    report, class_file = construct_filenames(kraken_outdir)
    class_file_data = read_kraken_file(class_file)
    report_info = read_kraken_report(report)
    report_data = associate_taxnames_with_readnames(report_info, class_file_data)
    write_custom_report(report_data, post_proc_report)
