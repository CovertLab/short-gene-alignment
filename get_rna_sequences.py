#! /usr/bin/env python
"""

"""

import argparse

def get_rna_sequences(rna_table_path, seq_path, output_path):
    with open(rna_table_path, 'r') as f:
        rna_ids = [line.split('\t')[0] for line in f.readlines()[1:]]

    rna_id_to_seq = {}
    with open(seq_path, 'r') as f:
        while len(rna_id_to_seq) < len(rna_ids):
            rna_id = f.readline()[:-1]
            if rna_id in rna_ids:
                seq = f.readline()[:-1]
                rna_id_to_seq[rna_id] = seq

    with open(output_path, 'w') as f:
        f.write('rna_id\tlength\tsequence\n')

        for rna_id, seq in rna_id_to_seq.items():
            f.write(f'{rna_id}\t{len(seq)}\t{seq}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract sequences of specific RNAs from a .seq file'
    )
    parser.add_argument(
        'rna_table_path', type=str,
        help='Path to the RNA table listing the IDs of RNAs to get sequences for'
    )
    parser.add_argument(
        'seq_path', type=str,
        help='Path to the .seq file listing the sequences of each RNA'
    )
    parser.add_argument(
        'output_path', type=str,
        help='Path and name of the output file'
    )

    args = parser.parse_args()

    get_rna_sequences(args.rna_table_path, args.seq_path, args.output_path)
