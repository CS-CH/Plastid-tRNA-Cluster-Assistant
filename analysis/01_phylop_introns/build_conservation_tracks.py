"""Build the all-species per-nucleotide phyloP/phastCons score track.
- Scores from 6 disjoint subsets are merged by alignment column: phyloP (CON/ACC) values
  present in every subset are averaged, and the per-subset values are retained.
- phastCons posterior probabilities are reported per reference position; they are mapped
  back to alignment columns through the non-gap columns of each subset's first row (row0).
- Coordinates: one reference sequence, chosen as the most complete across the intron loci.
Outputs:
  04_results/allspecies_track.npz
  04_results/perposition_<locus>.tsv
  04_results/region_summary_allspecies.tsv
  04_results/intron_occupancy.tsv
  04_results/reference_info.txt
"""
import os, csv, collections
import numpy as np

BASE = '/root/trna/work'
SRC = os.path.join(BASE, 'chloALL')
RES = '/mnt/c/Users/chang/Desktop/tRNA_migrate/04_results/task1_phylop_allspecies'
INTRON_LOCI = ['trnA-UGC', 'trnG-UCC', 'trnI-GAU', 'trnL-UAA']
SUBS = ['sub%d' % i for i in range(1, 7)]


def read_layout(p):
    out = []
    with open(p) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            out.append((r['locus'], int(r['col_start']), int(r['col_end'])))
    return out


def read_fasta(p):
    names, seqs = [], {}
    name, buf = None, []
    with open(p) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    names.append(name); seqs[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name:
        names.append(name); seqs[name] = ''.join(buf)
    return names, seqs


def read_wig(p):
    v = []
    with open(p) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('track') or s.startswith('fixedStep'):
                continue
            try:
                v.append(float(s))
            except ValueError:
                v.append(float('nan'))
    return np.array(v, dtype=float)


print('== Reading global alignment / layout / region labels ==')
gnames, gseqs = read_fasta(os.path.join(SRC, 'aln', 'concat.fasta'))
layout = read_layout(os.path.join(SRC, 'aln', 'layout.tsv'))
znp = np.load(os.path.join(SRC, 'col_region.npz'), allow_pickle=True)
region, cov = znp['region'], znp['cov']
NCOL = len(region)
print('  columns=%d  rows=%d' % (NCOL, len(gnames)))

# ---- Choose the reference: maximize the minimum non-gap count across the 4 intron loci ----
counts = collections.defaultdict(dict)
for locus, c0, c1 in layout:
    if locus not in INTRON_LOCI:
        continue
    for nm in gnames:
        s = gseqs[nm][c0 - 1:c1]
        counts[nm][locus] = sum(1 for ch in s if ch in 'ACGT')
best = None
for nm in gnames:
    d = counts[nm]
    if len(d) < len(INTRON_LOCI):
        continue
    score = (min(d.values()), sum(d.values()))
    if best is None or score > best[1]:
        best = (nm, score)
REF = best[0]
print('  reference =', REF, ' non-gap columns per locus =', counts[REF])

ref = gseqs[REF]

# ---- Merge scores from the 6 subsets ----
CON = np.full((len(SUBS), NCOL), np.nan)
ACC = np.full((len(SUBS), NCOL), np.nan)
PCsum = np.zeros(NCOL); PCn = np.zeros(NCOL)
for i, tag in enumerate(SUBS):
    W = os.path.join(BASE, tag)
    con = read_wig(os.path.join(W, 'out2', 'concat.CON.wig'))
    acc = read_wig(os.path.join(W, 'out2', 'concat.ACC.wig'))
    assert len(con) == NCOL, (tag, len(con), NCOL)
    CON[i] = con
    ACC[i] = acc
    pcf = sorted(f for f in os.listdir(os.path.join(W, 'pc2')) if f.startswith('pp_'))
    pc = read_wig(os.path.join(W, 'pc2', pcf[0]))
    names, seqs = read_fasta(os.path.join(W, 'aln', 'concat.fasta'))
    row0 = seqs[names[0]]
    # reference non-gap columns -> positions
    cols = np.array([j for j, ch in enumerate(row0) if ch != '-'])
    for k, j in enumerate(cols):
        if k < len(pc) and not np.isnan(pc[k]):
            PCsum[j] += pc[k]
            PCn[j] += 1
    print('  %s: CON/ACC ok; PC covers %d columns' % (tag, len(cols)))

mean_CON = np.nanmean(CON, axis=0)
mean_ACC = np.nanmean(ACC, axis=0)
with np.errstate(invalid='ignore', divide='ignore'):
    mean_PC = np.where(PCn > 0, PCsum / np.maximum(PCn, 1), np.nan)
np.savez(os.path.join(RES, 'allspecies_track.npz'),
         CON=CON, ACC=ACC, mean_CON=mean_CON, mean_ACC=mean_ACC,
         mean_PC=mean_PC, PCn=PCn, cov=cov,
         region=np.array([str(x) for x in region], dtype=object))
print('  saved track npz')


def slide(vals, w, min_frac=0.5):
    n = len(vals)
    out = np.full(n, np.nan)
    need = max(1, int(round(w * min_frac)))
    for i in range(n - w + 1):
        win = vals[i:i + w]
        v = win[~np.isnan(win)]
        if len(v) >= need:
            out[i + w // 2] = float(np.mean(v))
    return out


# ---- Per-locus per-nucleotide tables ----
os.makedirs(RES, exist_ok=True)
summary = []
for locus, c0, c1 in layout:
    if locus not in INTRON_LOCI:
        continue
    rs = ref[c0 - 1:c1]
    cols = np.array([j for j, ch in enumerate(rs) if ch in 'ACGT']) + (c0 - 1)
    v = mean_CON[cols]
    w15, w20, w30 = slide(v, 15), slide(v, 20), slide(v, 30)
    rows = []
    for k, j in enumerate(cols):
        rows.append(dict(
            locus=locus, ref_accession=REF, ref_pos=k + 1,
            align_col=int(j) + 1,
            ref_base=ref[j],
            region=str(region[j]),
            n_species_with_base=int(round(cov[j] * len(gnames))),
            coverage_3699='%.3f' % cov[j],
            phyloP_CON_mean='%.4f' % mean_CON[j] if not np.isnan(mean_CON[j]) else '',
            phyloP_ACC_mean='%.4f' % mean_ACC[j] if not np.isnan(mean_ACC[j]) else '',
            phastCons_mean='%.4f' % mean_PC[j] if not np.isnan(mean_PC[j]) else '',
            n_subsets_PC=int(PCn[j]),
            win15='%.4f' % w15[k] if not np.isnan(w15[k]) else '',
            win20='%.4f' % w20[k] if not np.isnan(w20[k]) else '',
            win30='%.4f' % w30[k] if not np.isnan(w30[k]) else '',
            **{'CON_' + t: ('%.4f' % CON[i][j] if not np.isnan(CON[i][j]) else '')
               for i, t in enumerate(SUBS)}))
    dst = os.path.join(RES, 'perposition_%s.tsv' % locus)
    with open(dst, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    print('  %s: %d reference positions -> %s' % (locus, len(rows), os.path.basename(dst)))
    # region summary
    for rg in ['intron', 'exon', 'flank']:
        m = np.array([str(region[j]) == rg for j in cols])
        idx = np.where(m)[0]
        if idx.size == 0:
            continue
        vals = v[idx]
        vals = vals[~np.isnan(vals)]
        if len(vals) == 0:
            continue
        summary.append(dict(locus=locus, region=rg, n_pos=int(idx.size),
                            CON_mean='%.4f' % np.mean(vals),
                            CON_median='%.4f' % np.median(vals),
                            frac_p_lt_0_05='%.3f' % float(np.mean(vals >= 1.30103)),
                            frac_p_lt_0_01='%.3f' % float(np.mean(vals >= 2.0)),
                            mean_coverage='%.3f' % np.mean([cov[j] for j in cols[idx]])))

with open(os.path.join(RES, 'region_summary_allspecies.tsv'), 'w', newline='', encoding='utf-8') as f:
    wr = csv.DictWriter(f, delimiter='\t', fieldnames=list(summary[0].keys()))
    wr.writeheader()
    wr.writerows(summary)

# ---- Intron occupancy ----
regdef = collections.defaultdict(dict)
with open(os.path.join(SRC, 'regions.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        regdef[r['locus']][r['acc']] = (int(r['gene_len']), int(r['intron_start']), int(r['intron_end']))
lines = ['locus\tn_in_alignment\tn_with_intron\tpct']
for locus in INTRON_LOCI:
    d = regdef[locus]
    tot = len(d)
    ni = sum(1 for a, v in d.items() if v[1] >= 0)
    lines.append('%s\t%d\t%d\t%.1f' % (locus, tot, ni, 100.0 * ni / max(1, tot)))
open(os.path.join(RES, 'intron_occupancy.tsv'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\nIntron occupancy:')
print('\n'.join(lines))

open(os.path.join(RES, 'reference_info.txt'), 'w', encoding='utf-8').write(
    'reference accession = %s\ncoverage per locus = %s\n' % (REF, counts[REF]))
print('\nReference sequence:', REF)
print('done')
