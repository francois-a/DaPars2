#!/usr/bin/python3
import os
import sys
import getopt
import tempfile
from collections import defaultdict
import argparse
import subprocess


def annotation_prepare_3UTR_extraction(gene_bed_file, gene_symbol_map_kfXref_file, output_utr_file):
    """"""
    with open(output_utr_file, 'w') as output_write:

        refseq_trapt_gene_symbol_dict = {}
        num_line = 0
        for line in open(gene_symbol_map_kfXref_file, 'r'):
            if num_line > 0:
                fields = line.strip('\n').strip('\r').split('\t')
                gene_symbol = fields[1]
                refseq_transcript_id = fields[0]
                refseq_trapt_gene_symbol_dict[refseq_transcript_id] = gene_symbol
            else:
                num_line += 1

        scanned_3UTR_list = []
        num_saved = 0
        for line in open(gene_bed_file, 'r'):
            fields = line.strip('\n').split('\t')
            refseq_id = fields[3]
            if '_' not in fields[0]:

                if refseq_id not in refseq_trapt_gene_symbol_dict:
                    gene_symbol = "NA"
                else:
                    gene_symbol = refseq_trapt_gene_symbol_dict[refseq_id]

                UTR_id = [refseq_id, gene_symbol, fields[0], fields[5]]
                UTR_id_new = '|'.join(UTR_id)
                curr_strand = fields[5]
                if curr_strand == "+":
                    UTR_end = fields[2]
                    gene_start = int(fields[1])
                    UTR_start = str(gene_start + int(fields[-1].strip(',').split(',')[-1]) + 1)  # 1base
                elif curr_strand == "-":
                    gene_start = int(fields[1])
                    UTR_start = str(gene_start + 1)  # 1base
                    UTR_end   = str(gene_start + int(fields[10].split(',')[0]))  # 1base, included

                this_UTR = fields[0] + UTR_start + UTR_end + curr_strand
                if this_UTR not in scanned_3UTR_list:
                    write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, '0', curr_strand]
                    output_write.writelines('\t'.join(write_line) + '\n')
                    scanned_3UTR_list.append(this_UTR)
                    num_saved += 1

    print(f"Total extracted 3' UTR: {num_saved}")


def UTRs_subtract_refine(UTRs_all):
    """"""
    strand_info = UTRs_all[0].strip('\n').split('\t')[-1]
    if strand_info == '+':
        all_pos = []
        for curr_line in UTRs_all:
            left_pos = curr_line.strip('\n').split('\t')[1]
            all_pos.append(int(left_pos))
        selected_UTR = UTRs_all[all_pos.index(min(all_pos))]
    else:
        all_pos = []
        for curr_line in UTRs_all:
            left_pos = curr_line.strip('\n').split('\t')[2]
            all_pos.append(int(left_pos))
        selected_UTR = UTRs_all[all_pos.index(max(all_pos))]
    return selected_UTR


def subtract_different_strand_overlap(input_gene_bed_file, output_utr_file):
    """"""
    with tempfile.NamedTemporaryFile() as f:
        # bedtools
        subprocess.check_call(f'subtractBed -a {input_gene_bed_file} -b {input_gene_bed_file} -S > {f.name}', shell=True)

        read_subtract_result_dict = defaultdict(list)
        for line in open(f.name, 'r'):
            transcript_id = line.split('\t')[3].split('|')[0]
            read_subtract_result_dict[transcript_id].append(line)

        output_utr_write = open(output_utr_file, 'w')
        for curr_trans_id in read_subtract_result_dict:
            curr_3UTRs = read_subtract_result_dict[curr_trans_id]
            if len(curr_3UTRs) == 1:
                output_utr_write.writelines(curr_3UTRs[0])
            else:
                output_utr_write.writelines(UTRs_subtract_refine(curr_3UTRs))
        output_utr_write.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DaPars: extract UTRs from annotation')
    parser.add_argument('-b', '--bed', required=True, help='Gene annotation in BED format')
    parser.add_argument('-s', '--symbol', required=True, help='Gene symbol annotation file')
    parser.add_argument('-o', '--ofile', required=True, help='Output file')
    args = parser.parse_args()

    print("Generating regions ...")
    with tempfile.NamedTemporaryFile() as f:
        annotation_prepare_3UTR_extraction(args.bed, args.symbol, f.name)
        subtract_different_strand_overlap(f.name, args.ofile)
    print("Finished")
