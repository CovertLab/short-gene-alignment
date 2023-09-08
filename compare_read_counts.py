#! /usr/bin/env python
"""
Compare read counts calculated from manual alignments to those predicted
from the manuscript.
"""

import argparse
from matplotlib import pyplot as plt
import numpy as np
import os


def compare_read_counts(rna_table_path):
    rna_id_to_actual_read_counts = {}

    for filename in os.listdir('results'):
        results_file = os.path.join('results', filename)

        if os.path.isfile(results_file) and results_file.endswith('.tsv'):
            with open(results_file, 'r') as f:
                lines = f.readlines()[1:]
                rna_ids = [line.split('\t')[0] for line in lines]
                read_counts = [float(line.split('\t')[3][:-1]) for line in lines]

                for rna_id, read_count in zip(rna_ids, read_counts):
                    if rna_id in rna_id_to_actual_read_counts:
                        rna_id_to_actual_read_counts[rna_id].append(read_count)
                    else:
                        rna_id_to_actual_read_counts[rna_id] = [read_count]

    for k, v in rna_id_to_actual_read_counts.items():
        rna_id_to_actual_read_counts[k] = np.mean(v)

    with open(rna_table_path, 'r') as f:
        lines = f.readlines()[1:]
        rna_id_to_predicted_read_counts = {
            line.split('\t')[0]: float(line.split('\t')[2][:-1]) for line in lines
        }

    rna_ids = [k for k in rna_id_to_actual_read_counts.keys()]
    actual_read_counts = np.array([
        rna_id_to_actual_read_counts[rna_id] for rna_id in rna_ids])
    predicted_read_counts = np.array([
        rna_id_to_predicted_read_counts[rna_id] for rna_id in rna_ids])

    fig = plt.figure(figsize=(6.35, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(
        np.zeros_like(predicted_read_counts),
        np.log10(actual_read_counts + 1),
        label='before adjustment', clip_on=False, s=12)
    ax.scatter(
        np.log10(predicted_read_counts + 1),
        np.log10(actual_read_counts + 1),
        label='after adjustment', clip_on=False, s=12)
    ax.set_xlim([0, 3.0])
    ax.set_ylim([0, 1.6])
    ax.set_xlabel('log10(TPM + 1)')
    ax.set_ylabel('log10(Manual alignment read counts + 1)')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 15))
    ax.spines["left"].set_position(("outward", 15))
    ax.legend(bbox_to_anchor=(1.1, 1))

    plt.tight_layout()
    plt.savefig('compare_read_counts.png')

    with open('rna_read_counts.tsv', 'w') as f:
        f.write('rna_id\tread_counts\n')
        for k, v in rna_id_to_actual_read_counts.items():
            f.write(f'{k}\t{v}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare read counts calculated from manual alignments to predicted read counts.'
    )
    parser.add_argument(
        'rna_table_path', type=str,
        help='Path to the RNA table listing the predicted read counts'
    )

    args = parser.parse_args()

    compare_read_counts(args.rna_table_path)
