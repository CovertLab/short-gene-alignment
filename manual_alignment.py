#! /usr/bin/env python
"""
Manually align each RNASeq read from a fasta file to short RNA sequences
"""

import argparse
from Bio import Align
import numpy as np
import os

SEQ_STARTERS = {'N', 'A', 'T', 'G', 'C'}


def manual_alignment(read_file_path, rna_file_path):
    read_file_basename = os.path.basename(read_file_path).split('.')[0]

    with open(rna_file_path, 'r') as f:
        lines = f.readlines()[1:]
        rna_ids = [line.split('\t')[0] for line in lines]
        rna_seqs = [line.split('\t')[2][:-1] for line in lines]

    rna_id_to_read_counts = {rna_id: 0 for rna_id in rna_ids}
    rna_id_to_effective_length = {rna_id: 0 for rna_id in rna_ids}
    read_lengths = []

    aligner = Align.PairwiseAligner()
    aligner.gap_score = -10
    aligner.end_gap_score = 0
    aligner.mismatch_score = -10

    with open(read_file_path, 'r') as f:
        line = f.readline()
        while line:
            if line[0] in SEQ_STARTERS:
                read_seq = line[:-1]
                read_seq_len = len(read_seq)
                read_lengths.append(read_seq_len)

                for rna_id, rna_seq in zip(rna_ids, rna_seqs):
                    rna_seq_len = len(rna_seq)
                    alignments = aligner.align(read_seq, rna_seq)

                    if alignments[0].score == min(read_seq_len, rna_seq_len):
                        rna_id_to_read_counts[rna_id] += 1

            line = f.readline()

    read_lengths = np.array(read_lengths)

    for rna_id, rna_seq in zip(rna_ids, rna_seqs):
        rna_length = len(rna_seq)
        effective_length = (np.abs(rna_length - read_lengths) + 1).mean()
        rna_id_to_effective_length[rna_id] = effective_length

    with open(os.path.join('results', f'{read_file_basename}_short_rna_read_counts.tsv'), 'w') as f:
        f.write('rna_id\tread_counts\teffective_length\tnormalized_read_counts\n')
        for rna_id, counts in rna_id_to_read_counts.items():
            el = rna_id_to_effective_length[rna_id]
            norm_counts = counts / el
            f.write(f'{rna_id}\t{counts}\t{el}\t{norm_counts}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Manually align short RNA sequences against RNA-Seq reads and get counts for perfect alignments'
    )
    parser.add_argument(
        'read_file_path', type=str,
        help='Path to the .fastq file containing the raw RNA-Seq reads'
    )
    parser.add_argument(
        'rna_file_path', type=str,
        help='Path to the file containing sequences of the short RNAs'
    )

    args = parser.parse_args()

    manual_alignment(args.read_file_path, args.rna_file_path)
