"""Extract intergenic lengths for three tRNA pairs and draw their distributions."""
import csv
import os
import re

import numpy as np
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

XL = 'data/Dataset5_phylogeny.xlsx'  # XLSX: pair sheets with accession, species and intergenic sequence
RES = 'results'  # length TSV and distribution plot in PDF
SHEETS = [
    ('trnG-trnR', 'G-R'),
    ('trnL-trnF', 'L-F'),
    ('trnI-trnA', 'I-A'),
]
ORDER = [
    ('G-R', 'trnG-UCC $-$ trnR-UCU', '#1b9e77'),
    ('L-F', 'trnL-UAA $-$ trnF-GAA', '#7570b3'),
    ('I-A', 'trnI-GAU $-$ trnA-UGC', '#d95f02'),
]


def clean_seq(value):
    return re.sub(r'[^ACGTN]', '', str(value or '').upper().replace('U', 'T'))


def read_lengths(workbook):
    rows = []
    for sheet, pair in SHEETS:
        ws = workbook[sheet]
        header = header_row = None
        for i, values in enumerate(ws.iter_rows(values_only=True, max_row=10), 1):
            if values and any(str(cell or '').strip().lower() == 'species' for cell in values):
                header, header_row = values, i
                break
        if header is None:
            raise ValueError(f'No header found in {sheet}')

        columns = {}
        for i, value in enumerate(header):
            label = str(value or '').strip().lower()
            if label == 'species':
                columns['species'] = i
            elif label.startswith('genbank'):
                columns['accession'] = i
            elif 'intergenic' in label and 'length' not in label and 'len' not in label:
                columns['intergenic'] = i
        if len(columns) != 3:
            raise ValueError(f'Missing accession, species or intergenic column in {sheet}')

        for values in ws.iter_rows(min_row=header_row + 1, values_only=True):
            species = values[columns['species']]
            accession = values[columns['accession']]
            sequence = clean_seq(values[columns['intergenic']])
            if species and accession and sequence:
                rows.append([pair, str(accession).strip(), str(species).strip(), len(sequence)])
    return rows


def plot_figure(input_tsv, out_dir=RES):
    os.makedirs(out_dir, exist_ok=True)
    lengths = {}
    with open(input_tsv, encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            lengths.setdefault(row['Pair'], []).append(int(row['Length']))

    bins = np.arange(0, 1010, 20)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for pair, label, color in ORDER:
        values = np.array(lengths[pair])
        ax.hist(values, bins=bins, histtype='stepfilled', alpha=0.35,
                color=color, edgecolor=color, lw=1.2)
        ax.hist(values, bins=bins, histtype='step', color=color, lw=1.2)
        median = float(np.median(values))
        ax.axvline(median, color=color, ls='--', lw=1.2,
                   label=f'{label}: n={len(values)}, median={median:.0f} bp')
    ax.set_xlabel('Length of the intergenic region between the two tRNAs (bp)')
    ax.set_ylabel('Number of species records')
    ax.set_xlim(0, 1000)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=10, loc='upper right')
    ax.set_title('Length distribution of the intergenic regions for the three tRNA pairs '
                 '(Dataset 5)', fontsize=12)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/intergenic_length_distribution.pdf')
    plt.close(fig)


def main():
    os.makedirs(RES, exist_ok=True)
    workbook = openpyxl.load_workbook(XL, read_only=True)
    rows = read_lengths(workbook)
    workbook.close()

    output = f'{RES}/intergenic_lengths_3pairs.tsv'
    with open(output, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['Pair', 'Accession', 'Species', 'Length'])
        writer.writerows(rows)
    plot_figure(output)


if __name__ == '__main__':
    main()
