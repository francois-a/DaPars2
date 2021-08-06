#!/usr/bin/python3
import pandas as pd
import numpy as np
import os
import sys
import getopt
import tempfile
from collections import defaultdict
import argparse
import subprocess


def annotation_prepare_3UTR_extraction(gene_bed_file, gene_symbol_map_kfXref_file, output_utr_file):
    """Parse UTRs (first/last exon depending on strand) from annotation"""

    # load gene symbol lookup
    refseq_trapt_gene_symbol_dict = pd.read_csv(gene_symbol_map_kfXref_file,
                sep='\t', index_col=0, squeeze=True).to_dict()

    with open(output_utr_file, 'w') as output_write:
        scanned_3UTR_list = []
        for line in open(gene_bed_file, 'r'):
            fields = line.strip('\n').split('\t')
            refseq_id = fields[3]
            if '_' not in fields[0]:  # skip alternate contigs

                gene_symbol = refseq_trapt_gene_symbol_dict.get(refseq_id, 'NA')
                # UTR ID: gene ID, gene name, chr, strand
                UTR_id = '|'.join([refseq_id, gene_symbol, fields[0], fields[5]])
                curr_strand = fields[5]
                gene_start = int(fields[1])
                # UTR: last exon of the gene
                if curr_strand == "+":
                    UTR_start = str(gene_start + int(fields[-1].strip(',').split(',')[-1]) + 1)  # 1base
                    UTR_end = fields[2]  # end of the gene
                elif curr_strand == "-":
                    UTR_start = str(gene_start + 1)  # 1base
                    UTR_end   = str(gene_start + int(fields[10].split(',')[0]))  # 1base, included

                # UTR: chr_start_end_strand --> only retain first occurrence
                this_UTR = fields[0] + UTR_start + UTR_end + curr_strand
                if this_UTR not in scanned_3UTR_list:
                    # chr, start, end, id/name/chr/strand, '0', strand
                    output_write.writelines('\t'.join([fields[0], UTR_start, UTR_end, UTR_id, '0', curr_strand]) + '\n')
                    scanned_3UTR_list.append(this_UTR)

    print(f"Total extracted 3' UTR: {len(scanned_3UTR_list)}")


def UTRs_subtract_refine(UTRs_all):
    """"""
    strand_info = UTRs_all[0].strip('\n').split('\t')[-1]
    if strand_info == '+':
        # first start position
        all_pos = np.array([int(curr_line.strip('\n').split('\t')[1]) for curr_line in UTRs_all])
        selected_UTR = UTRs_all[all_pos.argmin()]
    else:
        # last end position
        all_pos = np.array([int(curr_line.strip('\n').split('\t')[2]) for curr_line in UTRs_all])
        selected_UTR = UTRs_all[all_pos.argmax()]
    return selected_UTR


def subtract_different_strand_overlap(input_gene_bed_file, output_utr_file):
    """"""
    with tempfile.NamedTemporaryFile() as f:
        # subtract regions that overlap on opposite strands
        subprocess.check_call(f'bedtools subtract -a {input_gene_bed_file} -b {input_gene_bed_file} -S > {f.name}', shell=True)

        # group UTRs by refseq ID -- not optimal; overlapping transcripts in example data have different IDs
        read_subtract_result_dict = defaultdict(list)
        for line in open(f.name, 'r'):
            transcript_id = line.split('\t')[3].split('|')[0]  # UTR_id -> refseq_id
            read_subtract_result_dict[transcript_id].append(line)

        with open(output_utr_file, 'w') as output_utr_write:
            for _, curr_3UTRs in read_subtract_result_dict.items():
                if len(curr_3UTRs) == 1:
                    output_utr_write.writelines(curr_3UTRs[0])
                else:
                    output_utr_write.writelines(UTRs_subtract_refine(curr_3UTRs))


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
