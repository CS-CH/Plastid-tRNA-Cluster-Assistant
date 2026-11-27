"""Compare intron and flanking conservation across independent species subsets."""
import os, csv, math, random, collections
import numpy as np
from scipy import stats

RES = 'results'
INTRON_LOCI = ['trnA-UGC', 'trnG-UCC', 'trnI-GAU', 'trnL-UAA']
MINCOV = 0.5
NFRAME = 10
BLOCK = 30


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
            v.append(float(s))
    return np.array(v, dtype=float)


def ttest_1samp(x):
    x = np.asarray(x, dtype=float)
    n = len(x)
    m = x.mean()
    s = x.std(ddof=1)
    if s == 0 or n < 2:
        return m, float('nan'), float('nan')
    t = m / (s / math.sqrt(n))
    p = float(2 * stats.t.sf(abs(t), n - 1))
    return m, t, p


rng = random.Random(2026)
data = collections.defaultdict(lambda: collections.defaultdict(dict))

for i in range(1, 7):
    tag = 'sub%d' % i
    D = tag
    names, seqs = read_fasta(f'{D}/aln/concat.fasta')
    layout = read_layout(f'{D}/aln/layout.tsv')
    reg = np.load(f'{D}/col_region.npz', allow_pickle=True)['region']
    con = read_wig(f'{D}/out2/concat.CON.wig')
    M = len(names)
    for locus, c0, c1 in layout:
        if locus not in INTRON_LOCI:
            continue
        rows = [seqs[nm][c0 - 1:c1] for nm in names]
        ncols = c1 - c0 + 1
        occ = np.zeros(ncols)
        for s in rows:
            occ += np.isin(np.frombuffer(s.encode(), dtype=np.uint8),
                           np.frombuffer(b'ACGT', dtype=np.uint8))
        occ /= M
        reg_l, con_l = reg[c0 - 1:c1], con[c0 - 1:c1]
        dblock, dsite, mi_site, mf_site = [], [], [], []
        for nm in rng.sample(names, NFRAME):
            rs = seqs[nm][c0 - 1:c1]
            has = np.frombuffer(rs.encode(), dtype=np.uint8) != ord('-')
            idx = np.where((occ >= MINCOV) & has)[0]
            if len(idx) < 200:
                continue
            vv, rr = con_l[idx], reg_l[idx]
            # block-based contrast
            nb = len(vv) // BLOCK
            bi, bf = [], []
            for b in range(nb):
                sl = slice(b * BLOCK, (b + 1) * BLOCK)
                vals, rb = vv[sl], rr[sl]
                cnt = collections.Counter(rb)
                cnt.pop('gap', None)
                if not cnt:
                    continue
                lab, c = cnt.most_common(1)[0]
                if c < BLOCK * 0.5:
                    continue
                m = np.nanmean(vals)
                if np.isnan(m):
                    continue
                if lab == 'intron':
                    bi.append(m)
                elif lab == 'flank':
                    bf.append(m)
            if len(bi) >= 5 and len(bf) >= 5:
                dblock.append(float(np.mean(bi) - np.mean(bf)))
            # per-nucleotide contrast
            mi, mf = rr == 'intron', rr == 'flank'
            if mi.sum() >= 20 and mf.sum() >= 20:
                dsite.append(float(np.nanmean(vv[mi]) - np.nanmean(vv[mf])))
                mi_site.append(float(np.nanmean(vv[mi]))); mf_site.append(float(np.nanmean(vv[mf])))
        if not dblock:
            continue
        data[locus][tag] = dict(d_block=float(np.median(dblock)),
                                d_site=float(np.median(dsite)) if dsite else float('nan'),
                                int_mean=float(np.mean(mi_site)) if mi_site else float('nan'),
                                flk_mean=float(np.mean(mf_site)) if mf_site else float('nan'))

print('=== Hierarchical test: 6 disjoint species subsets = 6 independent estimates ===')
print('Within each subset: (intron mean - flank mean) over core columns (>=50% coverage),')
print('median over 10 reference frames')
print()
print('%-9s %-10s %10s %10s %8s %10s %10s' %
      ('locus', 'mode', 'intron_mean', 'flank_mean', 'diff', 't', 'p(2-tail)'))
rows_out = []
for locus in INTRON_LOCI:
    d = data[locus]
    if not d:
        continue
    tags = sorted(d)
    for mode in ['block', 'site']:
        vals = [d[t]['d_' + mode] for t in tags]
        mi_v = np.mean([d[t]['int_mean'] for t in tags])
        mf_v = np.mean([d[t]['flk_mean'] for t in tags])
        m, t, p = ttest_1samp(vals)
        print('%-9s %-10s %10.3f %10.3f %8.3f %10.2f %10.2e  subset diffs: %s' % (
            locus, 'block30bp' if mode == 'block' else 'per-site',
            mi_v, mf_v, m, t, p, ' '.join('%.2f' % v for v in vals)))
        rows_out.append(dict(locus=locus, mode=mode, mean_diff=m, t=t, p=p,
                             n_set=len(tags), all_positive=all(v > 0 for v in vals),
                             int_mean=mi_v, flk_mean=mf_v))
    print('-' * 100)

os.makedirs(RES, exist_ok=True)
dst = f'{RES}/phyloP_hierarchical_test.tsv'
with open(dst, 'w', newline='', encoding='utf-8') as f:
    wr = csv.DictWriter(f, delimiter='\t',
                        fieldnames=['locus', 'mode', 'n_set', 'int_mean', 'flk_mean',
                                    'mean_diff', 't', 'p', 'all_positive'])
    wr.writeheader()
    for r in rows_out:
        wr.writerow({k: ('%.4f' % v if isinstance(v, float) else v) for k, v in r.items()})
print()
print('saved', dst)
