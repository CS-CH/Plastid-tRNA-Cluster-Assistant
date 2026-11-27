"""Usage: python make_species_subsets.py <1..6> [K] [seed] [fit columns]."""
import os, sys, re, random, shutil

SRC = 'chloALL'  # taxa.txt, FASTA alignment, TSV layout, NPZ labels and Newick tree
gi = int(sys.argv[1])
K = int(sys.argv[2]) if len(sys.argv) > 2 else 6
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20261027
NFIT = int(sys.argv[4]) if len(sys.argv) > 4 else 1500


def parse_newick(s):
    s = s.strip()
    pos = 0

    def parse_node():
        nonlocal pos
        children = []
        if s[pos] == '(':
            pos += 1
            while True:
                children.append(parse_node())
                if s[pos] == ',':
                    pos += 1
                elif s[pos] == ')':
                    pos += 1
                    break
                else:
                    raise ValueError('unexpected %r at %d' % (s[pos], pos))
        m = re.match(r"[^,():;\[\]]*", s[pos:])
        name = m.group(0)
        pos += len(name)
        length = None
        if pos < len(s) and s[pos] == ':':
            pos += 1
            m = re.match(r"[-+0-9.eE]*", s[pos:])
            length = float(m.group(0))
            pos += len(m.group(0))
        return dict(children=children, name=name, length=length)

    return parse_node()


def write_node(n, out):
    if n['children']:
        out.append('(')
        for i, c in enumerate(n['children']):
            if i:
                out.append(',')
            write_node(c, out)
        out.append(')')
    out.append(n['name'] or '')
    if n['length'] is not None:
        out.append(':%.6f' % n['length'])


def tips(n, acc=None):
    acc = [] if acc is None else acc
    if not n['children']:
        acc.append(n['name'])
    for c in n['children']:
        tips(c, acc)
    return acc


def prune(n, keep):
    if not n['children']:
        return dict(children=[], name=n['name'], length=n['length']) if n['name'] in keep else None
    kids = [k for k in (prune(c, keep) for c in n['children']) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        c = kids[0]
        if c['length'] is not None and n['length'] is not None:
            c['length'] += n['length']
        return c
    return dict(children=kids, name=n['name'], length=n['length'])


tag = 'sub%d' % gi
dst = tag
for d in ['aln', 'tree', 'out2', 'log2']:
    os.makedirs(f'{dst}/{d}', exist_ok=True)

taxa = [l.strip() for l in open(f'{SRC}/taxa.txt') if l.strip()]
rng = random.Random(SEED)
order = taxa[:]
rng.shuffle(order)
group = set(a for i, a in enumerate(order) if i % K == gi - 1)
print('[%s] K=%d seed=%d -> %d species' % (tag, K, SEED, len(group)))

# ---- Alignment: keep subset rows, retain all columns ----
names, seqs, name, buf = [], {}, None, []
with open(f'{SRC}/aln/concat.fasta') as f:
    for line in f:
        if line.startswith('>'):
            if name:
                names.append(name); seqs[name] = ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
if name:
    names.append(name); seqs[name] = ''.join(buf)

keep = [n for n in names if n in group]
if len(keep) != len(group):
    print('  !! alignment rows (%d) != subset size (%d) -> missing %s' %
          (len(keep), len(group), sorted(group - set(names))[:3]))
ncol = len(seqs[names[0]])
print('  rows=%d cols=%d' % (len(keep), ncol))

with open(f'{dst}/aln/concat.fasta', 'w') as o:
    for n in keep:
        s = seqs[n]
        o.write('>' + n + '\n')
        for j in range(0, len(s), 80):
            o.write(s[j:j + 80] + '\n')
with open(f'{dst}/taxa.txt', 'w') as o:
    o.write('\n'.join(keep) + '\n')

good = []
for i in range(ncol):
    nb = 0
    for n in keep:
        if seqs[n][i] in 'ACGT':
            nb += 1
    if nb >= 0.5 * len(keep):
        good.append(i)
step = max(1, len(good) // NFIT)
sel = good[::step][:NFIT]
print('  fit cols: usable=%d -> keep=%d' % (len(good), len(sel)))
with open(f'{dst}/aln/fit.fasta', 'w') as o:
    for n in keep:
        s = seqs[n]
        o.write('>' + n + '\n' + ''.join(s[i] for i in sel) + '\n')

# ---- Tree pruning ----
root = parse_newick(open(f'{SRC}/tree/concat.treefile').read())
use = set(keep) & set(tips(root))
print('  tree target tips=%d' % len(use))
r = prune(root, use)
while r['children'] and len(r['children']) == 1:
    c = r['children'][0]
    if c['length'] is not None and r['length'] is not None:
        c['length'] += r['length']
    r = c
r['length'] = None
out = []
write_node(r, out)
nwk = ''.join(out)
open(f'{dst}/tree/pruned.treefile', 'w').write(nwk + '\n')
bang = re.sub(r':([0-9eE\.\-]+)', r':\1!', nwk)
open(f'{dst}/tree/bang.treefile', 'w').write(bang + '\n')
bl = [float(x) for x in re.findall(r':([0-9eE\.\-]+)!', bang)]
print('  pruned tips=%d branches=%d maxlen=%.4f' % (len(tips(r)), len(bl), max(bl)))

shutil.copy(f'{SRC}/aln/layout.tsv', f'{dst}/aln/layout.tsv')
shutil.copy(f'{SRC}/col_region.npz', f'{dst}/col_region.npz')
print('  done ->', dst)
