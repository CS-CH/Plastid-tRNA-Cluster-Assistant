"""Estimate intron-minus-flank effects and test them across species subsets."""
import os, csv, math
import numpy as np
from scipy import stats

RES = 'results'
INTRON_LOCI = ['trnA-UGC', 'trnG-UCC', 'trnI-GAU', 'trnL-UAA']
SUBS = ['sub%d' % i for i in range(1, 7)]
MINCOV = 0.5


def num(v):
    return np.nan if v == '' else float(v)


rows_out = []
stats_rows = []
for locus in INTRON_LOCI:
    rows = []
    with open(f'{RES}/perposition_{locus}.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            rows.append(r)
    diffs = []
    for t in SUBS:
        key = 'CON_' + t
        ci = [num(r[key]) for r in rows if r['region'] == 'intron' and num(r['coverage_3699']) >= MINCOV]
        cf = [num(r[key]) for r in rows if r['region'] == 'flank' and num(r['coverage_3699']) >= MINCOV]
        ce = [num(r[key]) for r in rows if r['region'] == 'exon' and num(r['coverage_3699']) >= MINCOV]
        ci = np.array([v for v in ci if not np.isnan(v)])
        cf = np.array([v for v in cf if not np.isnan(v)])
        ce = np.array([v for v in ce if not np.isnan(v)])
        d = float(np.mean(ci) - np.mean(cf)) if len(ci) and len(cf) else float('nan')
        diffs.append(d)
        rows_out.append(dict(locus=locus, subset=t,
                             n_intron=len(ci), n_flank=len(cf),
                             mean_intron='%.4f' % np.mean(ci) if len(ci) else '',
                             mean_flank='%.4f' % np.mean(cf) if len(cf) else '',
                             mean_exon='%.4f' % np.mean(ce) if len(ce) else '',
                             diff_intron_minus_flank='%.4f' % d if not np.isnan(d) else ''))
    diffs = np.array([d for d in diffs if not np.isnan(d)])
    m = float(np.mean(diffs))
    s = float(np.std(diffs, ddof=1))
    tstat = m / (s / math.sqrt(len(diffs)))
    p = float(2 * stats.t.sf(abs(tstat), len(diffs) - 1))
    stats_rows.append(dict(locus=locus, n_sets=len(diffs),
                           mean_diff='%.4f' % m, sd='%.4f' % s,
                           t='%.3f' % tstat, p='%.3e' % p,
                           all_six_positive='yes' if all(d > 0 for d in diffs) else 'no'))

with open(f'{RES}/subset_effects.tsv', 'w', newline='', encoding='utf-8') as f:
    wr = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows_out[0].keys()))
    wr.writeheader()
    wr.writerows(rows_out)
print('saved subset_effects.tsv')
print('%-9s %6s %10s %8s %8s %10s %6s' % ('locus', 'n', 'mean_diff', 'sd', 't', 'p', 'all+'))
for r in stats_rows:
    print('%-9s %6d %10s %8s %8s %10s %6s' % (r['locus'], r['n_sets'], r['mean_diff'],
                                              r['sd'], r['t'], r['p'], r['all_six_positive']))
with open(f'{RES}/hierarchical_test_from_effects.tsv', 'w', newline='', encoding='utf-8') as f:
    wr = csv.DictWriter(f, delimiter='\t', fieldnames=list(stats_rows[0].keys()))
    wr.writeheader()
    wr.writerows(stats_rows)
print('saved hierarchical_test_from_effects.tsv')
