"""Sequence features, conserved motifs and similarity of three tRNA intergenic regions."""
import os, re, math, csv, shutil, inspect, subprocess, collections
import numpy as np
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from openpyxl.styles import Font, Alignment, PatternFill

XL = 'data/Dataset5_phylogeny.xlsx'
RECORDS = 'data/records.tsv'
OUTROOT = 'results'
TABLES = f'{OUTROOT}/tables'
RAWDATA = f'{OUTROOT}/rawdata'
FIGURES = f'{OUTROOT}/figures'
CODE = f'{OUTROOT}/code'
XLSX = f'{OUTROOT}/Dataset_IG_analysis.xlsx'
TMP = f'{OUTROOT}/tmp'
OUT_TXT = f'{OUTROOT}/task2_intergenic_3pairs.txt'

K = 8
MIN_COV = 0.10          # present in >= 10% of records
ENT_MIN = 1.00          # low complexity: Shannon entropy < 1.0 bit
MAXRUN = 6              # low complexity: longest homopolymer >= 6 nt
PIQR_MAX = 0.15         # dispersed position: relative-position IQR > 0.15
ENR_MED = 20.0          # low-enrichment threshold
ENR_HIGH = 50.0         # high-reliability threshold
LEAK_PCT = 50.0         # possible flank-gene leakage threshold
NTOP = 12               # motifs shown per pair in position plots
NBIN = 20               # bins in position histograms
NT = 8

MIN_OV = 4              # minimum overlap for core-element assembly
KEEP_FRAC = 0.85        # minimum coverage relative to the smaller member
POS_TOL = 0.10          # maximum relative-position difference between merged 8-mers
MAXLEN = 60             # maximum core-element length
TOPN = 5                # entries retained per pair in the high-confidence shortlist
DEDUP = 0.05            # minimum relative-position difference between selected entries

# Three pairs only; G-fM excluded
SHEETS = [
    ('trnG-trnR', 'G-R', 'trnG-UCC', 'trnR-UCU'),
    ('trnI-trnA', 'I-A', 'trnI-GAU', 'trnA-UGC'),
    ('trnL-trnF', 'L-F', 'trnL-UAA', 'trnF-GAA'),
]
PAIRS = ['G-R', 'I-A', 'L-F']

TOB_RE = re.compile(r'nicotiana\s+tabacum', re.I)
MISSING = {'', '-', 'NAN', 'NONE', 'NA', '.', 'NF'}
TIER_COLOR = {'High': '#2b6cb0', 'Medium': '#dd8a2b', 'Low': '#c53030'}
TIER_ORDER = ['High', 'Medium', 'Low']


# ---------------------------------------------------------------- utilities
def norm_acc(s):
    return str(s or '').strip().split(',')[0].split(';')[0].split('.')[0].upper()


def clean_seq(s):
    if s is None:
        return ''
    t = str(s).strip().upper().replace('U', 'T').replace(' ', '')
    t = re.sub(r'[^ACGTN]', '', t)
    return '' if t in MISSING else t


_TR = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(_TR)[::-1]


def entropy(s):
    c = collections.Counter(s)
    n = len(s)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def maxrun(s):
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


def load_group():
    g = {}
    if not os.path.exists(RECORDS):
        return g
    with open(RECORDS, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            g[norm_acc(r['Accession'])] = r.get('Group', '') or ''
    return g


def read_sheet(wb, sheet):
    """Return (rows, ref_seq), where rows contain (acc, species, seq, flank)."""
    ws = wb[sheet]
    hdr, hdr_i = None, 0
    for i, row in enumerate(ws.iter_rows(values_only=True, max_row=10), 1):
        if row and any(str(c or '').strip().lower() == 'species' for c in row):
            hdr, hdr_i = row, i
            break
    if hdr is None:
        return [], ''
    ci = {}
    for j, h in enumerate(hdr):
        hl = str(h or '').strip().lower()
        if hl == 'species':
            ci['sp'] = j
        if hl.startswith('genbank'):
            ci['gb'] = j
        if 'intergenic' in hl and 'length' not in hl and 'len' not in hl:
            ci['ig'] = j
        if 'intron' in hl or 'spacer' in hl or hl.startswith('spliced'):
            ci.setdefault('flank', []).append(j)
    if any(k not in ci for k in ('sp', 'gb', 'ig')):
        return [], ''
    out, ref = [], ''
    for row in ws.iter_rows(min_row=hdr_i + 1, values_only=True):
        if not row or not row[ci['sp']]:
            continue
        sp = str(row[ci['sp']]).strip()
        acc = norm_acc(row[ci['gb']])
        seq = clean_seq(row[ci['ig']])
        if not acc or not seq:
            continue
        flank = ''.join(clean_seq(row[j]) for j in ci.get('flank', [])
                        if j < len(row))
        if TOB_RE.search(sp):
            ref = seq
        out.append((acc, sp, seq, flank))
    return out, ref


def reliability(lowc, piqr, enr, leak, cov_pct):
    """Assign a reliability tier and reason while retaining all motifs."""
    flags = []
    if lowc:
        flags.append('low-complexity repeat')
    if piqr > PIQR_MAX:
        flags.append(f'dispersed position (IQR={piqr:.2f})')
    if leak >= LEAK_PCT:
        flags.append(f'possible flank-gene leakage ({leak:.0f}% co-occurrence)')
    if enr < ENR_MED:
        flags.append(f'low enrichment ({enr:.0f}x)')
    if cov_pct < 15:
        flags.append(f'low coverage ({cov_pct:.0f}%)')

    if leak >= LEAK_PCT or (lowc and piqr > PIQR_MAX):
        tier = 'Low'
    elif lowc or piqr > PIQR_MAX:
        tier = 'Medium'
    elif enr >= ENR_HIGH:
        tier = 'High'
    else:
        tier = 'Medium'
    return tier, ('; '.join(flags) if flags else 'no anomaly')


def scan(seqs, pat):
    """Return plus-strand coverage, occurrences, positions and reverse-only records."""
    m = len(pat)
    r = rc(pat)
    cov = occ = rconly = 0
    pos = []
    for s in seqs:
        k = s.count(pat)
        if k:
            cov += 1
            occ += k
            den = len(s) - m
            if den <= 0:
                den = 1
            st = 0
            for _ in range(k):
                i = s.index(pat, st)
                pos.append(i / den)
                st = i + 1
        elif r in s:
            rconly += 1
    return cov, occ, pos, rconly


def try_merge(a, b):
    """Return the largest-overlap merge of a and b, including containment."""
    for x in range(min(len(a), len(b)), MIN_OV - 1, -1):
        j = b.find(a[-x:])                     # suffix of a lands at position j of b
        if j < 0:
            continue
        s = j + x - len(a)                     # start of a in b's coordinate frame
        lo, hi = max(0, s), min(len(b), s + len(a))
        if a[lo - s:hi - s] != b[lo:hi]:       # overlapping bases must match exactly
            continue
        return (b[:s] if s > 0 else '') + a + (b[hi:] if hi < len(b) else '')
    for x in range(min(len(a), len(b)), MIN_OV - 1, -1):
        j = a.find(b[-x:])                     # reverse orientation
        if j < 0:
            continue
        s = j + x - len(b)
        lo, hi = max(0, s), min(len(a), s + len(b))
        if b[lo - s:hi - s] != a[lo:hi]:
            continue
        return (a[:s] if s > 0 else '') + b + (a[hi:] if hi < len(a) else '')
    return None


def assemble(rows, seqs):
    """Assemble overlapping High 8-mers and verify each merge against real sequences."""
    units = [dict(contig=r['km'], members=[r]) for r in rows]
    changed = True
    while changed:
        changed = False
        for i in range(len(units)):
            for j in range(len(units)):
                if i == j:
                    continue
                ui, uj = units[i], units[j]
                if abs(float(ui['members'][0]['pmed'])
                       - float(uj['members'][0]['pmed'])) > POS_TOL:
                    continue
                new = try_merge(ui['contig'], uj['contig'])
                if new is None or len(new) > MAXLEN:
                    continue
                cov = scan(seqs, new)[0]
                if cov < KEEP_FRAC * min(int(ui['members'][0]['cov']),
                                         int(uj['members'][0]['cov'])):
                    continue
                units[i] = dict(contig=new,
                                members=ui['members'] + uj['members'])
                units.pop(j)
                changed = True
                break
            if changed:
                break
    return units


# ------------------------------------------- figure data loaders (standalone)
def load_fig1(raw):
    """Read intergenic length values per pair for Figure 1."""
    d = collections.defaultdict(list)
    with open(f'{raw}/fig1c_length_values.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            d[r['Pair']].append(int(r['Length_bp']))
    return {k: np.array(v) for k, v in d.items()}


def load_fig2(raw):
    """Read length and GC values per pair for Figure 2."""
    d = collections.defaultdict(lambda: ([], []))
    with open(f'{raw}/fig2b_length_vs_gc.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            x, y = d[r['Pair']]
            x.append(int(r['Length_bp']))
            y.append(float(r['GC_pct']))
    return {k: (np.array(v[0]), np.array(v[1])) for k, v in d.items()}


def load_fig3(raw):
    """Read conserved 8-mer position summaries for Figure 3."""
    d = collections.defaultdict(list)
    with open(f'{raw}/fig3c_motif_summary.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            d[r['Pair']].append(dict(
                km=r['Kmer'], tier=r['Reliability'],
                pmed=float(r['Pos_median']),
                tob=float(r['TobaccoCoord_bp']),
                tob10=float(r['Tobacco_p10']), tob25=float(r['Tobacco_p25']),
                tob75=float(r['Tobacco_p75']), tob90=float(r['Tobacco_p90'])))
    return d


def load_fig4(raw):
    """Read BLASTN HSP details per pair for Figure 4."""
    d = collections.defaultdict(list)
    with open(f'{raw}/fig4b_blast_hits.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            d[r['Pair']].append((r['Query'], r['Subject'],
                                 float(r['Pident']), int(r['Length'])))
    return d


# -------------------------------------------------------------- plotting
def plot_fig1(L, out, pairs):
    """Intergenic length distributions."""
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(5.0 * n, 4.2))
    axes = np.atleast_1d(axes)
    for ax, pr in zip(axes, pairs):
        v = L[pr]
        ax.hist(v, bins=40, color='steelblue', edgecolor='black', lw=.4)
        ax.axvline(np.median(v), color='crimson', ls='--', lw=1.5,
                   label=f'median = {np.median(v):.0f} bp')
        ax.set_xlabel('Intergenic length (bp)')
        ax.set_ylabel('Records')
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_fig2(gc, tiers, out, pairs):
    """GC boxplots, length versus GC, and motif reliability tiers."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    bp = axes[0].boxplot([gc[pr][1] for pr in pairs],
                         showfliers=False, patch_artist=True)
    axes[0].set_xticks(range(1, len(pairs) + 1))
    axes[0].set_xticklabels(pairs)
    for b in bp['boxes']:
        b.set_facecolor('#8fbcd4')
    axes[0].set_ylabel('GC content (%)')
    axes[0].set_title('(a) GC content by pair')
    for pr in pairs:
        axes[1].scatter(gc[pr][0], gc[pr][1], s=6, alpha=.25, label=pr)
    axes[1].set_xlabel('Intergenic length (bp)')
    axes[1].set_ylabel('GC content (%)')
    axes[1].set_title('(b) Length vs GC')
    axes[1].legend(fontsize=8, markerscale=2)
    cnt = collections.Counter()
    for pr in pairs:
        for d in tiers.get(pr, []):
            cnt[d['tier']] += 1
    axes[2].bar(TIER_ORDER, [cnt.get(t, 0) for t in TIER_ORDER],
                color=[TIER_COLOR[t] for t in TIER_ORDER])
    axes[2].set_ylabel('Number of 8-mers')
    axes[2].set_title(f'(c) Reliability of conserved 8-mers (top-{NTOP}/pair)')
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_fig3(tiers, ref_len, out, pairs):
    """Positions of conserved 8-mers on the tobacco reference scale."""
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(5.6 * n, 5.2))
    axes = np.atleast_1d(axes)
    for ax, pr in zip(axes, pairs):
        recs = sorted(tiers.get(pr, []), key=lambda d: d['pmed'])
        for i, d in enumerate(recs):
            c = TIER_COLOR[d['tier']]
            ax.hlines(i, d['tob10'], d['tob90'], color=c, lw=1, alpha=.7)
            ax.hlines(i, d['tob25'], d['tob75'], color=c, lw=7, alpha=.9)
            ax.scatter([d['tob']], [i], color=c, s=30, zorder=3)
        ax.set_yticks(range(len(recs)))
        ax.set_yticklabels([f'{d["km"]}  [{d["tier"]}]' for d in recs],
                           fontsize=8, family='monospace')
        ax.set_ylim(-1, len(recs))
        Lref = ref_len[pr]
        ax.axvline(0, color='k', ls=':', lw=1)
        ax.axvline(Lref, color='k', ls=':', lw=1)
        ax.set_xlabel('Position in intergenic region (bp; tobacco reference scale)')
        ax.set_title(f'{pr}  (tobacco spacer = {Lref} bp)', fontsize=10)
    handles = [Patch(color=TIER_COLOR[t], label=f'{t} confidence')
               for t in TIER_ORDER]
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=10)
    fig.suptitle('Task 2 - Position of conserved 8-mers within the intergenic region\n'
                 '(0 = immediately downstream of upstream gene; L = immediately upstream of downstream gene)',
                 fontsize=13)
    fig.tight_layout(rect=(0, .05, 1, .93))
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_fig4(hits, out, pairs):
    """Cross-species identity and alignment-length distributions."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    for pr in pairs:
        hs = hits.get(pr, [])
        if hs:
            axes[0].hist([h[2] for h in hs], bins=np.arange(80, 100.5, 1),
                         alpha=.55, label=pr)
    axes[0].set_xlabel('BLASTN percent identity (%)')
    axes[0].set_ylabel('Number of HSPs')
    axes[0].set_title('(a) Cross-species identity of intergenic regions')
    axes[0].legend(fontsize=8)
    for pr in pairs:
        hs = hits.get(pr, [])
        if hs:
            axes[1].hist([h[3] for h in hs], bins=40, alpha=.55, label=pr)
    axes[1].set_xlabel('Alignment length (bp)')
    axes[1].set_ylabel('Number of HSPs')
    axes[1].set_title('(b) Alignment length distribution')
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


# ------------------------------------------------------------- code export
CODE_HEADER = '''# -*- coding: utf-8 -*-
"""{title}

Generated by: analyze_intergenic_regions.py (three-pair analysis; G-fM excluded)
Purpose: reproduce {png} independently from the rawdata tables
Usage: python {fname}
Requirements: numpy, matplotlib
"""
import os, csv, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RAW = r"{raw}"
FIG = r"{fig}"

TIER_COLOR = {{'High': '#2b6cb0', 'Medium': '#dd8a2b', 'Low': '#c53030'}}
TIER_ORDER = ['High', 'Medium', 'Low']
NTOP = {ntop}
'''

CODE_FOOTER = '''

if __name__ == '__main__':
    out = FIG + '/{png}'
    {plot}({args}, out, PAIRS)
    print('saved:', out)
'''


def export_code(codedir, raw_win, fig_win, pairs):
    """Export each plotting function as an independently runnable script."""
    os.makedirs(codedir, exist_ok=True)
    made = []

    specs = [
        (1, 'Figure 1: intergenic length distributions', 'fig1_length_distribution.png',
         'plot_fig1', "load_fig1(RAW)", ['load_fig1']),
        (2, 'Figure 2: GC composition and motif reliability', 'fig2_gc_and_reliability.png',
         'plot_fig2', "load_fig2(RAW), load_fig3(RAW)",
         ['load_fig2', 'load_fig3']),
        (3, 'Figure 3: positions of conserved 8-mers', 'fig3_motif_position.png',
         'plot_fig3', "load_fig3(RAW), REF_LEN", ['load_fig3']),
        (4, 'Figure 4: cross-species sequence similarity', 'fig4_similarity.png',
         'plot_fig4', "load_fig4(RAW)", ['load_fig4']),
    ]
    for num, title, png, plotname, argstr, loaders in specs:
        fname = f'make_fig{num}.py'
        parts = [CODE_HEADER.format(title=title, png=png, fname=fname,
                                    raw=raw_win, fig=fig_win, ntop=NTOP)]
        parts.append('PAIRS = %r' % (pairs,))
        if num == 3:
            parts.append('REF_LEN = %r' % (ref_len_store,))
        for fn in loaders:
            parts.append(inspect.getsource(globals()[fn]))
        parts.append(inspect.getsource(globals()[plotname]))
        parts.append(CODE_FOOTER.format(plot=plotname, args=argstr, png=png))
        p = f'{codedir}/{fname}'
        with open(p, 'w', encoding='utf-8') as f:
            f.write('\n\n'.join(parts))
        made.append(p)

    # Include a copy of the complete workflow for traceability
    dst = f'{codedir}/{os.path.basename(__file__)}'
    if os.path.abspath(dst) != os.path.abspath(__file__):
        shutil.copy2(__file__, dst)

    readme = f'{codedir}/README.md'
    with open(readme, 'w', encoding='utf-8') as f:
        f.write('\n'.join([
            '# Plotting code (three pairs: G-R / I-A / L-F)',
            '',
            'These scripts were generated by `analyze_intergenic_regions.py`.',
            'They use only `rawdata/*.tsv` and can be run independently of the main workflow.',
            '',
            '| Script | Output | Input data |',
            '|---|---|---|',
            '| `make_fig1.py` | `../figures/fig1_length_distribution.png` | `fig1c_length_values.tsv` |',
            '| `make_fig2.py` | `../figures/fig2_gc_and_reliability.png` | `fig2b_length_vs_gc.tsv`, `fig3c_motif_summary.tsv` |',
            '| `make_fig3.py` | `../figures/fig3_motif_position.png` | `fig3c_motif_summary.tsv` |',
            '| `make_fig4.py` | `../figures/fig4_similarity.png` | `fig4b_blast_hits.tsv` |',
            '',
            '```bash',
            'python make_fig1.py   # equivalent for the other figures',
            '```',
            '',
            '`analyze_intergenic_regions.py` is the complete workflow that generated these scripts.',
            '',
            '**Scope**: only G-R, I-A and L-F are included;',
            'G-fM is excluded because it is a monocot-only subset without a tobacco reference.',
            '',
        ]))
    return made


ref_len_store = {}
blast_store = {}          # pair -> [(q, s, pident, length)]


# ------------------------------------------------------------------- Excel
def make_excel(data, tiers, sp_rows, comp_rows, mot_rows, sim_rows,
               high_rows, rd_len, rd_gc, rd_blast, rd_len_gc, rd_possum,
               rd_blast_raw, spg, xlsx, BINS, GCBINS, PIDBINS):
    from openpyxl.chart import BarChart, ScatterChart, Reference, Series
    from openpyxl.chart.marker import Marker

    pairs = [p for p in PAIRS if p in data]
    wb = openpyxl.Workbook()

    def put(title, header, rows, freeze=True):
        ws = wb.create_sheet(title)
        ws.append(header)
        for c in ws[1]:
            c.font = Font(bold=True)
            c.fill = PatternFill('solid', fgColor='DDEBF7')
            c.alignment = Alignment(vertical='center')
        for r in rows:
            ws.append(r)
        for col in ws.columns:
            w = max((len(str(c.value)) for c in col if c.value is not None),
                    default=0)
            ws.column_dimensions[col[0].column_letter].width = \
            min(max(9, w + 2), 45)
        if freeze:
            ws.freeze_panes = 'A2'
        return ws

    # ---- README ----
    ws = wb.active
    ws.title = 'README'
    rows = [
        ['Task 2: intergenic-region analysis of three tRNA pairs (G-fM excluded)'],
        [''],
        ['Data source',
         'data/Dataset5_phylogeny.xlsx ("Dataset 1" in the reviewer comment)'],
        ['Extraction',
         'Read the "Intergenic trnX-XXX/trnY-YYY" column from each pair sheet'],
        ['Pairs',
         'G-R (trnG-UCC/trnR-UCU), I-A (trnI-GAU/trnA-UGC), '
         'L-F (trnL-UAA/trnF-GAA)'],
        ['Reason for excluding G-fM',
         'G-fM (trnG-UCC/trnfM-CAU) is a monocot-only subset (n=354), much smaller than the other three; '
         'it lacks Nicotiana tabacum, preventing use of a common tobacco coordinate scale; '
         'it is therefore not comparable with the other three pairs and is excluded'],
        [''],
        ['Length (spacing)',
         'Intergenic lengths are summarized and grouped as <50 / 50-100 / 100-200 / 200-300 / '
         '300-500 / 500-1k / >=1k bp'],
        ['Composition',
         'Per-record GC%, total GC and A/C/G/T frequencies; AT-skew = (A-T)/(A+T)'],
        ['Conserved motifs',
         f'k={K}; record coverage (N_records/Pct_records); '
         'Enrichment = observed occurrences / occurrences expected under an independent mononucleotide model; '
         'Position = relative coordinate from 0 (next to upstream gene) to 1 (next to downstream gene)'],
        ['Tobacco reference coordinates',
         'The Nicotiana tabacum (AP019624) intergenic length provides the reference scale; '
         'relative positions are mapped linearly to bp (all three pairs include tobacco)'],
        ['Sequence similarity',
         'All-against-all BLASTN self-search (evalue<=1e-5, dust off, word_size 7, len>=20bp)'],
        [''],
        ['=== Motif reliability tiers (all motifs retained) ==='],
        ['High (reliable)', 'No anomaly: not low complexity + concentrated position (IQR<=0.15) + enrichment >=50x'],
        ['Medium (interpret with caution)',
         'Low-complexity repeat, enrichment of 20-50x, or coverage <15%'],
        ['Low (not recommended for interpretation)',
         'Possible flank-gene leakage (>=50% co-occurrence in spliced tRNA/intron columns), '
         'or dispersed position (IQR>0.15) combined with low complexity'],
        ['Note column', 'Reports the values triggering each classification for record-level review'],
        [''],
        ['=== High-confidence shortlist ==='],
        ['Method',
         'Use only High 8-mers; greedily assemble sequence overlaps (>=%d nt) into core elements; '
         'verify every merge against real sequences and reject it if coverage falls below %.0f%% of the smaller member; '
         'sort by coverage and enrichment, then deduplicate by position (|Delta RelPos|<=%.2f); '
         'retain the top %d entries per pair' % (MIN_OV, KEEP_FRAC * 100, DEDUP, TOPN)],
        ['Enrichment definitions',
         'Enrichment_8mer_x uses the k=%d model and is directly comparable with T4_Motifs; '
         'Enrichment_element_x describes the complete core element and is for internal reference only'
         % K],
        [''],
        ['=== Worksheet navigation ==='],
        ['T1_Spacing', 'Intergenic length statistics and bin counts per pair'],
        ['T2_SpacingByGroup', 'Length statistics grouped by taxonomy'],
        ['T3_Composition', 'Base composition per pair (GC / A / C / G / T)'],
        ['T4_Motifs', f'All conserved {K}-mers with reliability tiers and notes'],
        ['T4b_HighConf', 'High-confidence core-element shortlist'],
        ['T5_Similarity', 'BLASTN cross-species similarity summary'],
        ['RD_Fig1', 'Source data for Figure 1 length histograms'],
        ['RD_Fig2', 'Source data for Figure 2 length versus GC scatter plots'],
        ['RD_Fig3', 'Source data for Figure 3 conserved 8-mer positions'],
        ['RD_Fig4', 'Source data for Figure 4 BLASTN identity and alignment length'],
        ['figures/', 'PNG figures produced with matplotlib'],
        ['code/', 'Standalone plotting code (make_fig1~4.py)'],
        [''],
    ]
    for r in rows:
        ws.append(r)
    ws['A1'].font = Font(bold=True, size=13)
    for r in range(3, ws.max_row + 1):
        ws.cell(row=r, column=1).font = Font(bold=True)
    ws.column_dimensions['A'].width = 22
    ws.column_dimensions['B'].width = 115
    for row in ws.iter_rows(min_row=3, min_col=2, max_col=2):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical='top')

    # Result tables
    put('T1_Spacing',
        ['Pair', 'Upstream', 'Downstream', 'n', 'min', '25%', 'median', '75%',
         'max', 'mean', 'sd'] + [l for _, _, l in BINS], sp_rows)
    put('T2_SpacingByGroup',
        ['Pair', 'Group', 'n', 'min', '25%', 'median', '75%', 'max', 'mean', 'sd'],
        [[k[0], k[1]] + list(v) for k, v in sorted(spg.items())])
    put('T3_Composition',
        ['Pair', 'Upstream', 'Downstream', 'n', 'len_min', 'len_median',
         'len_mean', 'len_max', 'GC_total%', 'GC_median%', 'A%', 'C%', 'G%', 'T%'],
        comp_rows)
    put('T4_Motifs',
        ['Pair', 'Upstream', 'Downstream', 'Kmer', 'N_records', 'Pct_records',
         'Occ_obs', 'Occ_exp', 'Enrichment_OccRatio', 'Pos_median', 'Pos_p10',
         'Pos_p25', 'Pos_p75', 'Pos_p90', 'Pos_IQR', 'Entropy_bits', 'MaxRun_nt',
         'LowComplexity', 'LeakPct', 'Reliability', 'Note', 'TobaccoCoord_bp'],
        mot_rows)
    put('T4b_HighConf',
        ['Pair', 'Upstream', 'Downstream', 'Rank', 'Core_element_5to3',
         'Len_nt', 'Core_8mer', 'N_records', 'N_total', 'Pct_records',
         'Pct_records_either_strand', 'Occ_obs', 'Occ_exp',
         'Enrichment_8mer_x', 'Enrichment_element_x', 'Pos_median',
         'Pos_IQR', 'TobaccoCoord_bp', 'Spacer_ref_bp', 'Entropy_bits',
         'MaxRun_nt', 'N_8mers_merged', 'Member_8mers', 'Reliability'],
        high_rows)
    put('T5_Similarity',
        ['Pair', 'Upstream', 'Downstream', 'N_seq', 'N_hits', 'N_seq_with_hit',
         'Pct_with_hit'], sim_rows)

    # Source data for figures
    lbl_bins = [l for _, _, l in BINS]
    piv = collections.defaultdict(dict)
    for pr, lbl, c, pct in rd_len:
        piv[lbl][pr] = c
    ws1 = put('RD_Fig1', ['Bin_bp'] + pairs,
              [[b] + [piv[b].get(p, 0) for p in pairs] for b in lbl_bins])
    ch = BarChart()
    ch.type = 'col'
    ch.grouping = 'clustered'
    ch.title = 'Intergenic length distribution (spacing)'
    ch.y_axis.title = 'Number of records'
    ch.x_axis.title = 'Intergenic length (bp)'
    ch.add_data(Reference(ws1, min_col=2, max_col=1 + len(pairs),
                          min_row=1, max_row=1 + len(lbl_bins)),
                titles_from_data=True)
    ch.set_categories(Reference(ws1, min_col=1, min_row=2,
                                max_row=1 + len(lbl_bins)))
    ch.width, ch.height = 20, 11
    ws1.add_chart(ch, 'H2')

    ws2 = put('RD_Fig2', ['Pair', 'Accession', 'Species', 'Group', 'Length_bp',
                          'GC_pct'], rd_len_gc)
    sc = ScatterChart()
    sc.title = 'Intergenic length vs GC content'
    sc.x_axis.title = 'Intergenic length (bp)'
    sc.y_axis.title = 'GC content (%)'
    sc.style = 13
    r0 = 2
    rng = np.random.default_rng(0)
    for pr in pairs:
        nrow = len([1 for _, _, s, _ in data[pr]['rows'] if s])
        idx = list(range(nrow))
        if len(idx) > 1200:
            idx = list(rng.choice(idx, 1200, replace=False))
        xs = Reference(ws2, min_col=5, min_row=r0, max_row=r0 + len(idx) - 1)
        ys = Reference(ws2, min_col=6, min_row=r0, max_row=r0 + len(idx) - 1)
        s = Series(ys, xs, title=pr)
        s.marker = Marker(symbol='circle', size=3)
        s.graphicalProperties.line.noFill = True
        sc.series.append(s)
        r0 += nrow
    sc.width, sc.height = 18, 11
    ws2.add_chart(sc, 'I2')

    ws3 = put('RD_Fig3',
              ['Pair', 'Kmer', 'Reliability', 'Pos_median', 'Pos_p10', 'Pos_p25',
               'Pos_p75', 'Pos_p90', 'TobaccoCoord_bp', 'Tobacco_p10',
               'Tobacco_p25', 'Tobacco_p75', 'Tobacco_p90'], rd_possum)
    for r in range(2, ws3.max_row + 1):
        t = ws3.cell(row=r, column=3).value
        if t in TIER_COLOR:
            ws3.cell(row=r, column=3).fill = PatternFill(
                'solid', fgColor=TIER_COLOR[t].lstrip('#'))
    row = 2
    for pr in pairs:
        nrow = len(tiers.get(pr, []))
        if not nrow:
            continue
        first = row
        row += nrow
        bch = BarChart()
        bch.type = 'bar'
        bch.title = (f'{pr}: conserved 8-mer position '
                     f'(tobacco spacer = {ref_len_store[pr]} bp)')
        bch.x_axis.title = 'Position in intergenic region (bp)'
        bch.add_data(Reference(ws3, min_col=9, min_row=first, max_row=row - 1),
                     titles_from_data=False)
        bch.set_categories(Reference(ws3, min_col=2, min_row=first,
                                     max_row=row - 1))
        bch.width, bch.height = 20, 12
        ws3.add_chart(bch, f'P{first}')

    lbl_p = [l for _, _, l in PIDBINS]
    piv4 = collections.defaultdict(dict)
    for pr, lbl, c, pct in rd_blast:
        piv4[lbl][pr] = c
    ws4b = put('RD_Fig4_pivot', ['Pident_bin_%'] + pairs,
               [[b] + [piv4[b].get(p, 0) for p in pairs] for b in lbl_p])
    ch4 = BarChart()
    ch4.type = 'col'
    ch4.grouping = 'clustered'
    ch4.title = 'BLASTN percent identity of intergenic regions (cross-species)'
    ch4.y_axis.title = 'Number of HSPs'
    ch4.x_axis.title = 'Percent identity (%)'
    ch4.add_data(Reference(ws4b, min_col=2, max_col=1 + len(pairs),
                           min_row=1, max_row=1 + len(lbl_p)),
                 titles_from_data=True)
    ch4.set_categories(Reference(ws4b, min_col=1, min_row=2,
                                 max_row=1 + len(lbl_p)))
    ch4.width, ch4.height = 20, 11
    ws4b.add_chart(ch4, 'H2')
    put('RD_Fig4', ['Pair', 'Query', 'Subject', 'Pident', 'Length'],
        rd_blast_raw)

    wb.save(xlsx)
    print('excel  ->', xlsx)


def blast_self(fa, tag):
    db = f'{TMP}/{tag}_db'
    subprocess.run(['makeblastdb', '-in', fa, '-dbtype', 'nucl', '-out', db],
                   check=True, capture_output=True)
    r = subprocess.run(['blastn', '-query', fa, '-db', db, '-evalue', '1e-5',
                        '-dust', 'no', '-word_size', '7', '-num_threads', str(NT),
                        '-outfmt', '6 qseqid sseqid pident length',
                        '-max_target_seqs', '20'], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError((r.stderr or 'blastn failed')[:300])
    hits = []
    for line in r.stdout.splitlines():
        p = line.split('\t')
        if len(p) < 4:
            continue
        q, s, pid, ln = p[0], p[1], float(p[2]), int(p[3])
        if q == s or ln < 20:
            continue
        hits.append((q.split('|')[0], s.split('|')[0], pid, ln))
    return hits


def main():
    for d in (TABLES, RAWDATA, FIGURES, CODE, TMP):
        os.makedirs(d, exist_ok=True)
    out = []
    P = out.append
    P('=== Task 2: three-pair analysis (G-R / I-A / L-F), G-fM excluded ===')
    P('  Length / composition / conserved motifs with reliability tiers / similarity / high-confidence shortlist')
    P('  Data source: data/Dataset5_phylogeny.xlsx')
    P('')

    grp = load_group()
    wb = openpyxl.load_workbook(XL, read_only=True)
    data = {}
    for sheet, pair, up, dn in SHEETS:
        rows, ref = read_sheet(wb, sheet)
        if rows:
            data[pair] = dict(up=up, dn=dn, rows=rows, ref=ref)
            P(f'  [{sheet}] n={len(rows)}  tobacco reference intergenic length='
              f'{len(ref) if ref else "N/A"} bp')
    wb.close()
    P('')

    BINS = [(0, 50, '<50'), (50, 100, '50-100'), (100, 200, '100-200'),
            (200, 300, '200-300'), (300, 500, '300-500'),
            (500, 1000, '500-1k'), (1000, 10 ** 9, '>=1k')]
    GCBINS = [(0, 15, '<15'), (15, 25, '15-25'), (25, 35, '25-35'),
              (35, 45, '35-45'), (45, 100, '>=45')]
    PIDBINS = [(80, 85, '80-85'), (85, 90, '85-90'), (90, 95, '90-95'),
               (95, 99, '95-99'), (99, 100.01, '>=99')]

    sp_rows, comp_rows, mot_rows, sim_rows, high_rows = [], [], [], [], []
    rd_len, rd_gc, rd_pos, rd_posbin, rd_blast = [], [], [], [], []
    rd_possum, rd_blast_raw = [], []
    spg = {}
    tiers = {}
    mot_all = {}          # pair -> all conserved 8-mers, not only the top N

    for pair in PAIRS:
        if pair not in data:
            continue
        up, dn = data[pair]['up'], data[pair]['dn']
        recs, ref = data[pair]['rows'], data[pair]['ref']
        seqs = [r[2] for r in recs]
        n = len(seqs)
        L_ref = len(ref) if ref else int(np.median([len(s) for s in seqs]))
        ref_len_store[pair] = L_ref
        P(f'########## {pair} ({up} / {dn}) ##########')
        P(f'  n={n}   tobacco reference length={L_ref} bp'
          f'{"" if ref else "  [no tobacco record; median length used]"}')

        # Length
        lens = np.array([len(s) for s in seqs])
        P(f'  Length (bp): min={lens.min()}  25%={np.percentile(lens,25):.0f}  '
          f'median={np.median(lens):.0f}  75%={np.percentile(lens,75):.0f}  '
          f'mean={lens.mean():.0f}  sd={lens.std():.1f}  max={lens.max()}')
        bcnt = []
        for a, b, lbl in BINS:
            c = int(((lens >= a) & (lens < b)).sum())
            bcnt.append(c)
            rd_len.append([pair, lbl, c, round(c / n * 100, 2)])
        for (a, b, lbl), c in zip(BINS, bcnt):
            P(f'    {lbl:>8} bp: {c:>6} ({c/n*100:.1f}%)')
        sp_rows.append([pair, up, dn, n, lens.min(), np.percentile(lens, 25),
                        np.median(lens), np.percentile(lens, 75), lens.max(),
                        round(float(lens.mean()), 1), round(float(lens.std()), 1)]
                       + bcnt)

        # Composition
        allseq = ''.join(seqs)
        base = {b: allseq.count(b) / len(allseq) * 100 for b in 'ACGT'}
        gc = base['G'] + base['C']
        at = base['A'] + base['T']
        gcv = np.array([(s.count('G') + s.count('C')) / len(s) * 100
                        for s in seqs if s])
        P(f'  GC%: 25%={np.percentile(gcv,25):.0f}  median={np.median(gcv):.0f}  '
          f'75%={np.percentile(gcv,75):.0f}  mean={gcv.mean():.1f}  '
          f'sd={gcv.std():.1f}')
        P('  Mononucleotides: ' + '  '.join(f'{b}={base[b]:.1f}%' for b in 'ACGT'))
        P(f'  Total GC={gc:.1f}%  AT-skew={(base["A"]-base["T"])/at:+.3f}')
        comp_rows.append([pair, up, dn, n, int(lens.min()),
                          int(np.median(lens)), round(float(lens.mean()), 1),
                          int(lens.max()), round(gc, 1),
                          round(float(np.median(gcv)), 1)]
                         + [round(base[b], 1) for b in 'ACGT'])
        for a, b, lbl in GCBINS:
            c = int(((gcv >= a) & (gcv < b)).sum())
            rd_gc.append([pair, lbl, c, round(c / len(gcv) * 100, 2)])

        # ---------- k-mer ----------
        occ = collections.Counter()
        cov = collections.Counter()
        pos = collections.defaultdict(list)
        for s in seqs:
            if len(s) < K:
                continue
            den = len(s) - K
            for i in range(len(s) - K + 1):
                km = s[i:i + K]
                if 'N' in km:
                    continue
                occ[km] += 1
                pos[km].append(i / den if den else 0.0)
            for km in {s[i:i + K] for i in range(len(s) - K + 1)
                       if 'N' not in s[i:i + K]}:
                cov[km] += 1

        p = {b: base[b] / 100 for b in 'ACGT'}
        tot_pos = sum(max(0, len(s) - K + 1) for s in seqs)
        keep = [km for km, c in cov.items() if c >= MIN_COV * n]

        rec = []
        for km in keep:
            e = tot_pos
            for ch in km:
                e *= p.get(ch, 0.25)
            enr = occ[km] / e if e > 0 else 0
            ra = np.array(pos[km])
            pmed = float(np.median(ra))
            p10 = float(np.percentile(ra, 10))
            p25 = float(np.percentile(ra, 25))
            p75 = float(np.percentile(ra, 75))
            p90 = float(np.percentile(ra, 90))
            piqr = p75 - p25
            ent, mr = entropy(km), maxrun(km)
            lowc = (ent < ENT_MIN) or (mr >= MAXRUN)
            iig = ifl = both = 0
            for _, _, s, fl in recs:
                a = km in s
                b = bool(fl) and km in fl
                iig += a
                ifl += b
                both += (a and b)
            leak = (both / iig * 100) if iig else 0
            tier, note = reliability(lowc, piqr, enr, leak, cov[km] / n * 100)
            tb = lambda x: round(x * L_ref, 1)
            rec.append(dict(km=km, cov=cov[km], pct=round(cov[km] / n * 100, 2),
                            occ=occ[km], exp=round(e, 2), enr=round(enr, 1),
                            pmed=round(pmed, 3), p10=round(p10, 3),
                            p25=round(p25, 3), p75=round(p75, 3),
                            p90=round(p90, 3), piqr=round(piqr, 3),
                            ent=round(ent, 2), mr=mr, lowc=lowc,
                            leak=round(leak, 1), tier=tier, note=note,
                            tob=tb(pmed), tob10=tb(p10), tob25=tb(p25),
                            tob75=tb(p75), tob90=tb(p90)))
        rec.sort(key=lambda d: (-d['cov'], -d['enr']))

        n_hi = sum(1 for d in rec if d['tier'] == 'High')
        P(f'  --- Conserved {K}-mers (coverage>={MIN_COV*100:.0f}%), {len(rec)} total; '
          f'all retained and tiered (High={n_hi}) ---')
        P(f'  {"kmer":<10}{"records":>6}{"%":>7}{"enrich":>8}{"pos":>7}{"IQR":>7}'
          f'{"entropy":>6}{"run":>5}  tier  note')
        for _i, d in enumerate(rec):
            if _i < NTOP:
                P(f'  {d["km"]:<10}{d["cov"]:>6}{d["pct"]:>7}{d["enr"]:>8}'
                  f'{d["pmed"]:>7}{d["piqr"]:>7}{d["ent"]:>6}{d["mr"]:>5}'
                  f'  {d["tier"]:<6}  {d["note"]}')
            mot_rows.append([pair, up, dn, d['km'], d['cov'], d['pct'],
                             d['occ'], d['exp'], d['enr'], d['pmed'],
                             d['p10'], d['p25'], d['p75'], d['p90'], d['piqr'],
                            d['ent'], d['mr'], 'yes' if d['lowc'] else 'no',
                            d['leak'], d['tier'], d['note'], d['tob']])

        top = rec[:NTOP]
        tier_of = {d['km']: d['tier'] for d in top}
        for acc, sp, s, _ in recs:
            L = len(s)
            if L < K:
                continue
            den = L - K
            for i in range(L - K + 1):
                km = s[i:i + K]
                if km in tier_of:
                    rel = i / den if den else 0.0
                    rd_pos.append([pair, km, acc, sp, grp.get(acc, ''), L, i,
                                   round(rel, 4), round(rel * L_ref, 1), i,
                                   L - i - K, tier_of[km]])
        bstep = 1.0 / NBIN
        for d in rec[:NTOP]:
            ra = np.array(pos[d['km']])
            cnt, _ = np.histogram(ra, bins=NBIN, range=(0, 1))
            for b in range(NBIN):
                if cnt[b]:
                    rd_posbin.append([pair, d['km'], d['tier'],
                                      round((b + 0.5) * bstep, 3),
                                      round((b + 0.5) * bstep * L_ref, 1),
                                      int(cnt[b])])
            rd_possum.append([pair, d['km'], d['tier'], d['pmed'], d['p10'],
                              d['p25'], d['p75'], d['p90'], d['tob'],
                              d['tob10'], d['tob25'], d['tob75'], d['tob90']])
        tiers[pair] = rec[:NTOP]
        mot_all[pair] = rec
        P('')

        # Similarity
        fa = f'{TMP}/{pair}.fasta'
        with open(fa, 'w', encoding='utf-8') as f:
            for i, (acc, sp, s, _) in enumerate(recs):
                f.write(f'>{acc}|{i}\n')
                for j in range(0, len(s), 70):
                    f.write(s[j:j + 70] + '\n')
        hits = blast_self(fa, pair)
        blast_store[pair] = hits
        with_sp = len({h[0] for h in hits})
        P('  --- BLASTN self-search (e<=1e-5, dust off, word_size 7, len>=20bp) ---')
        P(f'    Significant HSPs={len(hits)}  with cross-species match={with_sp}/{n} '
          f'({with_sp/n*100:.1f}%)')
        if hits:
            hl = np.array([h[3] for h in hits])
            hp = np.array([h[2] for h in hits])
            P(f'    Alignment length: median={np.median(hl):.0f}bp  max={hl.max()}bp')
            P(f'    Identity: median={np.median(hp):.1f}%  mean={hp.mean():.1f}%  '
              f'>=90% fraction={(hp>=90).mean()*100:.1f}%')
            for a, b, lbl in PIDBINS:
                c = int(((hp >= a) & (hp < b)).sum())
                rd_blast.append([pair, lbl, c, round(c / len(hp) * 100, 2)])
            for h in hits:
                rd_blast_raw.append([pair, h[0], h[1], h[2], h[3]])
        sim_rows.append([pair, up, dn, n, len(hits), with_sp,
                         round(with_sp / n * 100, 1)])
        P('')

    # Length by taxonomy group
    for pair, cfg in data.items():
        by = collections.defaultdict(list)
        for acc, sp, s, _ in cfg['rows']:
            by[grp.get(acc, '') or '(unknown)'].append(len(s))
        for g, v in by.items():
            a = np.array(v)
            spg[(pair, g)] = (len(a), a.min(), np.percentile(a, 25),
                              np.median(a), np.percentile(a, 75), a.max(),
                              round(float(a.mean()), 1), round(float(a.std()), 1))

    # High-confidence shortlist
    P('########## High-confidence shortlist: High 8-mers assembled into core elements ##########')
    for pair in PAIRS:
        if pair not in data:
            continue
        up, dn = data[pair]['up'], data[pair]['dn']
        seqs = [r[2] for r in data[pair]['rows']]
        n = len(seqs)
        L_ref = ref_len_store[pair]
        allseq = ''.join(seqs)
        p = {b: allseq.count(b) / len(allseq) for b in 'ACGT'}
        sub = sorted([d for d in mot_all[pair] if d['tier'] == 'High'],
                     key=lambda d: (-d['cov'], -d['enr']))
        units = assemble(sub, seqs) if sub else []
        cands = []
        for u in units:
            contig = u['contig']
            m = len(contig)
            cov, occ, posl, rconly = scan(seqs, contig)
            if not cov:
                continue
            tot = sum(max(0, len(s) - m + 1) for s in seqs)
            e = tot
            for ch in contig:
                e *= p.get(ch, 0.25)
            ra = np.array(posl)
            pmed = float(np.median(ra))
            p25, p75 = float(np.percentile(ra, 25)), float(np.percentile(ra, 75))
            best = min(u['members'], key=lambda d: (-d['cov'], -d['enr']))
            cands.append(dict(
                pair=pair, up=up, dn=dn, core=best['km'], elem=contig, m=m,
                cov=cov, n=n, pct=round(cov / n * 100, 2),
                pct_rc=round((cov + rconly) / n * 100, 2),
                occ=occ, exp=round(e, 4),
                enr_e=round(occ / e, 1) if e > 0 else 0,
                enr8=best['enr'], pmed=round(pmed, 3),
                piqr=round(p75 - p25, 3), ent=round(entropy(contig), 2),
                mr=maxrun(contig), nkm=len(u['members']),
                members=' '.join(x['km'] for x in u['members']),
                L_ref=L_ref, tob=round(pmed * L_ref, 1)))
        cands.sort(key=lambda d: (-d['cov'], -d['enr8']))
        top = []
        for c in cands:
            if any(try_merge(s['elem'], c['elem']) is not None
                   or abs(s['pmed'] - c['pmed']) <= DEDUP for s in top):
                continue
            top.append(c)
            if len(top) >= TOPN:
                break
        for i, t in enumerate(top, 1):
            t = dict(t, rank=i)
            high_rows.append([t['pair'], t['up'], t['dn'], t['rank'], t['elem'],
                              t['m'], t['core'], t['cov'], t['n'], t['pct'],
                              t['pct_rc'], t['occ'], t['exp'], t['enr8'],
                              t['enr_e'], t['pmed'], t['piqr'], t['tob'],
                              t['L_ref'], t['ent'], t['mr'], t['nkm'],
                              t['members'], 'High'])
        P(f'  {pair}: High 8-mer={len(sub)} -> core element={len(cands)} '
          f'-> selected {len(top)}')
    P('')

    # Write TSV tables
    def dump(sub, name, header, rows):
        p = (RAWDATA if sub == 'raw' else TABLES) + '/' + name
        with open(p, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(header)
            w.writerows(rows)
        return p

    dump('tab', 'ig_spacing.tsv',
         ['Pair', 'Upstream', 'Downstream', 'n', 'min', '25%', 'median', '75%',
          'max', 'mean', 'sd'] + [l for _, _, l in BINS], sp_rows)
    dump('tab', 'ig_composition.tsv',
         ['Pair', 'Upstream', 'Downstream', 'n', 'len_min', 'len_median',
          'len_mean', 'len_max', 'GC_total%', 'GC_median%', 'A%', 'C%', 'G%', 'T%'],
         comp_rows)
    dump('tab', 'ig_motifs_annotated.tsv',
         ['Pair', 'Upstream', 'Downstream', 'Kmer', 'N_records', 'Pct_records',
          'Occ_obs', 'Occ_exp', 'Enrichment_OccRatio', 'Pos_median', 'Pos_p10',
          'Pos_p25', 'Pos_p75', 'Pos_p90', 'Pos_IQR', 'Entropy_bits', 'MaxRun_nt',
          'LowComplexity', 'LeakPct', 'Reliability', 'Note',
          'TobaccoCoord_bp'], mot_rows)
    dump('tab', 'ig_similarity.tsv',
         ['Pair', 'Upstream', 'Downstream', 'N_seq', 'N_hits', 'N_seq_with_hit',
          'Pct_with_hit'], sim_rows)
    dump('tab', 'ig_spacing_by_group.tsv',
         ['Pair', 'Group', 'n', 'min', '25%', 'median', '75%', 'max', 'mean', 'sd'],
         [[k[0], k[1]] + list(v) for k, v in sorted(spg.items())])
    dump('tab', 'ig_motifs_highconf.tsv',
         ['Pair', 'Upstream', 'Downstream', 'Rank', 'Core_element_5to3',
          'Len_nt', 'Core_8mer', 'N_records', 'N_total', 'Pct_records',
          'Pct_records_either_strand', 'Occ_obs', 'Occ_exp',
          'Enrichment_8mer_x', 'Enrichment_element_x', 'Pos_median',
          'Pos_IQR', 'TobaccoCoord_bp', 'Spacer_ref_bp', 'Entropy_bits',
          'MaxRun_nt', 'N_8mers_merged', 'Member_8mers', 'Reliability'],
         high_rows)

    with open(f'{TABLES}/ig_details.tsv', 'w', newline='',
              encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Pair', 'Accession', 'Species', 'Group', 'Length', 'GC%',
                    'Sequence'])
        for pair, cfg in data.items():
            for acc, sp, s, _ in cfg['rows']:
                g = (s.count('G') + s.count('C')) / len(s) * 100 if s else 0
                w.writerow([pair, acc, sp, grp.get(acc, ''), len(s), f'{g:.1f}', s])

    dump('raw', 'fig1a_length_bins.tsv', ['Pair', 'Bin_bp', 'Count', 'Pct'], rd_len)
    dump('raw', 'fig1b_length_by_group.tsv',
         ['Pair', 'Group', 'n', 'min', '25%', 'median', '75%', 'max', 'mean', 'sd'],
         [[k[0], k[1]] + list(v) for k, v in sorted(spg.items())])
    dump('raw', 'fig1c_length_values.tsv',
         ['Pair', 'Accession', 'Species', 'Group', 'Length_bp'],
         [[pair, acc, sp, grp.get(acc, ''), len(s)]
          for pair, cfg in data.items() for acc, sp, s, _ in cfg['rows']])
    dump('raw', 'fig2a_gc_bins.tsv', ['Pair', 'GC_bin_%', 'Count', 'Pct'], rd_gc)
    dump('raw', 'fig2b_length_vs_gc.tsv',
         ['Pair', 'Accession', 'Species', 'Group', 'Length_bp', 'GC_pct'],
         [[pair, acc, sp, grp.get(acc, ''), len(s),
           round((s.count('G') + s.count('C')) / len(s) * 100, 2) if s else 0]
          for pair, cfg in data.items() for acc, sp, s, _ in cfg['rows']])
    dump('raw', 'fig3a_motif_occurrences.tsv',
         ['Pair', 'Kmer', 'Accession', 'Species', 'Group', 'SeqLen_bp',
          'Pos_bp', 'RelPos_0_1', 'TobaccoCoord_bp', 'DistToUpstream_bp',
          'DistToDownstream_bp', 'Reliability'], rd_pos)
    dump('raw', 'fig3b_motif_position_bins.tsv',
         ['Pair', 'Kmer', 'Reliability', 'RelPos_bin', 'TobaccoCoord_bp',
          'Count'], rd_posbin)
    dump('raw', 'fig3c_motif_summary.tsv',
         ['Pair', 'Kmer', 'Reliability', 'Pos_median', 'Pos_p10', 'Pos_p25',
          'Pos_p75', 'Pos_p90', 'TobaccoCoord_bp', 'Tobacco_p10',
          'Tobacco_p25', 'Tobacco_p75', 'Tobacco_p90'],
         [[pr, d['km'], d['tier'], d['pmed'], d['p10'], d['p25'], d['p75'],
           d['p90'], d['tob'], d['tob10'], d['tob25'], d['tob75'], d['tob90']]
          for pr in PAIRS if pr in tiers for d in tiers[pr]])
    dump('raw', 'fig4a_blast_pident_bins.tsv',
         ['Pair', 'Pident_bin_%', 'Count', 'Pct'], rd_blast)
    dump('raw', 'fig4b_blast_hits.tsv',
         ['Pair', 'Query', 'Subject', 'Pident', 'Length'], rd_blast_raw)

    # Figures
    pairs = [p for p in PAIRS if p in data]
    rd_len_gc = [[pr, acc, sp, grp.get(acc, ''), len(s),
                  round((s.count('G') + s.count('C')) / len(s) * 100, 2)
                  if s else 0]
                 for pr in pairs for acc, sp, s, _ in data[pr]['rows']]
    L = {pr: np.array([len(s) for _, _, s, _ in data[pr]['rows']])
         for pr in pairs}
    gcxy = {pr: (L[pr], np.array([(s.count('G') + s.count('C')) / len(s) * 100
                                 for _, _, s, _ in data[pr]['rows'] if s]))
            for pr in pairs}
    plot_fig1(L, f'{FIGURES}/fig1_length_distribution.png', pairs)
    plot_fig2(gcxy, tiers, f'{FIGURES}/fig2_gc_and_reliability.png',
              pairs)
    plot_fig3(tiers, ref_len_store, f'{FIGURES}/fig3_motif_position.png', pairs)
    plot_fig4(blast_store, f'{FIGURES}/fig4_similarity.png', pairs)
    print('figures ->', FIGURES)

    made = export_code(CODE, OUTROOT + r'\rawdata',
                       OUTROOT + r'\figures', pairs)
    print('code    ->', CODE, f'({len(made)} plotting scripts + a copy of the main script)')

    make_excel(data, tiers, sp_rows, comp_rows, mot_rows, sim_rows,
               high_rows, rd_len, rd_gc, rd_blast, rd_len_gc, rd_possum,
               rd_blast_raw, spg, XLSX, BINS, GCBINS, PIDBINS)

    txt = '\n'.join(out)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    with open(OUT_TXT, 'w', encoding='utf-8') as f:
        f.write(txt + '\n')
    print(txt)
    print('\nOutput directory:', OUTROOT)


if __name__ == '__main__':
    main()
