"""Combine subset phyloP CON scores and draw conservation profiles."""
import csv
import collections
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SRC = 'chloALL'  # FASTA alignment, TSV layout and NPZ coverage/region labels
RES = 'results'  # per-position TSV tables and conservation plot in PDF
LOCI = ['trnG-UCC', 'trnL-UAA', 'trnI-GAU', 'trnA-UGC']
SUBS = [f'sub{i}' for i in range(1, 7)]
COL = {'intron': '#d95f02', 'exon': '#1b9e77', 'flank': '#7570b3'}


def read_layout(path):
    with open(path) as handle:
        return [(r['locus'], int(r['col_start']), int(r['col_end']))
                for r in csv.DictReader(handle, delimiter='\t')]


def read_fasta(path):
    names, seqs = [], {}
    name, parts = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name:
                    names.append(name)
                    seqs[name] = ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name:
        names.append(name)
        seqs[name] = ''.join(parts)
    return names, seqs


def read_wig(path):
    values = []
    with open(path) as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith(('track', 'fixedStep')):
                values.append(float(value))
    return np.array(values, dtype=float)


def read_tsv(path):
    with open(path, encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def num(value):
    return np.nan if value is None or value == '' else float(value)


def slide(values, width=15):
    result = np.full(len(values), np.nan)
    required = max(1, int(round(width * 0.5)))
    for i in range(len(values) - width + 1):
        window = values[i:i + width]
        observed = window[~np.isnan(window)]
        if len(observed) >= required:
            result[i + width // 2] = np.mean(observed)
    return result


def spans(rows, region):
    result, start, previous = [], None, None
    for row in rows:
        position = int(row['ref_pos'])
        if row['region'] == region:
            if start is None:
                start = position
            previous = position
        elif start is not None:
            result.append((start, previous))
            start = None
    if start is not None:
        result.append((start, previous))
    return result


def merge(regions, gap=3):
    result = []
    for start, end in regions:
        if result and start - result[-1][1] - 1 <= gap:
            result[-1] = (result[-1][0], end)
        else:
            result.append((start, end))
    return result


def plot_figure(data_dir=RES, out_dir=RES):
    os.makedirs(out_dir, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8))
    for k, (ax, locus) in enumerate(zip(axes.ravel(), LOCI)):
        rows = read_tsv(f'{data_dir}/perposition_{locus}.tsv')
        x = np.array([int(r['ref_pos']) for r in rows])
        score = np.array([num(r['phyloP_CON_mean']) for r in rows])
        window = np.array([num(r['win15']) for r in rows])

        for region in ('flank', 'exon', 'intron'):
            for start, end in merge(spans(rows, region)):
                ax.axvspan(start, end, color=COL[region],
                           alpha=0.10 if region == 'flank' else 0.22, lw=0)
        ax.plot(x, score, color='#b0b0b0', lw=0.7)
        ax.plot(x, window, color='black', lw=1.5)

        exons = merge(spans(rows, 'exon'))
        labels = [('exon 1', exons[:1]),
                  ('intron', merge(spans(rows, 'intron'))),
                  ('exon 2', exons[1:] if len(exons) > 1 else [])]
        for label, regions in labels:
            for start, end in regions:
                if end - start >= 20:
                    ax.text((start + end) / 2, 1.015, label,
                            transform=ax.get_xaxis_transform(), ha='center',
                            va='bottom', fontsize=8, color='#333333')
        ax.set_title(locus, fontsize=12, pad=26)
        ax.set_xlabel('nucleotide position across the region (nt)')
        ax.set_ylabel('phyloP ($-\\log_{10}P$)')
        ax.grid(alpha=0.2, lw=0.5)
        ax.set_xlim(x.min(), x.max())
        ax.text(-0.09, 1.08, 'ABCD'[k], transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top')

    handles = [
        Patch(facecolor=COL['intron'], alpha=0.4, label='intron'),
        Patch(facecolor=COL['exon'], alpha=0.4, label='exon (tRNA)'),
        Patch(facecolor=COL['flank'], alpha=0.35, label='150-bp flanks'),
        plt.Line2D([], [], color='black', lw=1.5, label='15-bp sliding-window mean'),
        plt.Line2D([], [], color='#b0b0b0', lw=1.0, label='per-nucleotide phyloP'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=10, frameon=False)
    fig.suptitle('Per-nucleotide phyloP conservation across the four intron-containing '
                 'plastid tRNAs (mean of six species sets, 3,699 accessions)', fontsize=12.5)
    fig.tight_layout(rect=[0, 0.055, 1, 0.955], h_pad=2.2)
    fig.savefig(f'{out_dir}/phyloP_conservation_15bp_sliding_mean.pdf')
    plt.close(fig)


def main():
    os.makedirs(RES, exist_ok=True)
    names, sequences = read_fasta(f'{SRC}/aln/concat.fasta')
    layout = read_layout(f'{SRC}/aln/layout.tsv')
    labels = np.load(f'{SRC}/col_region.npz', allow_pickle=True)
    coverage, regions = labels['cov'], labels['region']

    counts = collections.defaultdict(dict)
    for locus, start, end in layout:
        if locus in LOCI:
            for name in names:
                seq = sequences[name][start - 1:end]
                counts[name][locus] = sum(base in 'ACGT' for base in seq)
    reference = max((name for name in names if len(counts[name]) == len(LOCI)),
                    key=lambda name: (min(counts[name].values()), sum(counts[name].values())))
    ref_seq = sequences[reference]

    subset_scores = []
    for subset in SUBS:
        score = read_wig(f'{subset}/out2/concat.CON.wig')
        assert len(score) == len(regions), (subset, len(score), len(regions))
        subset_scores.append(score)
    mean_score = np.nanmean(np.array(subset_scores), axis=0)

    for locus, start, end in layout:
        if locus not in LOCI:
            continue
        ref_part = ref_seq[start - 1:end]
        columns = np.array([i for i, base in enumerate(ref_part) if base in 'ACGT']) + start - 1
        score = mean_score[columns]
        window = slide(score)
        rows = []
        for i, column in enumerate(columns):
            rows.append({
                'locus': locus,
                'ref_accession': reference,
                'ref_pos': i + 1,
                'align_col': int(column) + 1,
                'ref_base': ref_seq[column],
                'region': str(regions[column]),
                'n_species_with_base': int(round(coverage[column] * len(names))),
                'coverage_3699': f'{coverage[column]:.3f}',
                'phyloP_CON_mean': '' if np.isnan(score[i]) else f'{score[i]:.4f}',
                'win15': '' if np.isnan(window[i]) else f'{window[i]:.4f}',
            })
        with open(f'{RES}/perposition_{locus}.tsv', 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, delimiter='\t', fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    plot_figure()


if __name__ == '__main__':
    main()
