import os
import sys
import json


def construct_filenames(kraken_outdir):
    class_file = os.path.join(kraken_outdir, "classifications.kraken")
    report = os.path.join(kraken_outdir, "classifications.report")
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
                kraken_class[class_status][seq_id] = info
    print(kraken_class)
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
    print(final_class_lines)
    return final_class_lines


def associate_taxnames_with_contigs(
    report_info: dict, class_file_data: dict, sort_key: str = "seq_len"  # "contig_id"
):
    report = []
    classified_tigs = class_file_data["C"]
    unclassified = class_file_data["U"]
    for tig in classified_tigs:
        tig_id = tig
        tig_info = classified_tigs[tig]
        tig_info["contig_id"] = tig_id
        for class_line in report_info:
            if class_line["ncbi-tax-id"] == tig_info["tax-id"]:
                tig_info["class_rank"] = class_line["taxon_rank"]
                tig_info["class_name"] = class_line["scientific_name"]
                report.append(tig_info)
    unclass_msg = f"{len(unclassified)} contigs could not be classified"
    report = sorted(report, key=lambda x: x[sort_key])[::-1]
    if len(unclassified):
        report.append(unclass_msg)
    print("report", report)
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
    report_data = associate_taxnames_with_contigs(report_info, class_file_data)
    write_custom_report(report_data, post_proc_report)
