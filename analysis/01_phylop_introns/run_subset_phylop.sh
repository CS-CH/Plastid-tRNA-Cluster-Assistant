#!/usr/bin/env bash
# Usage: bash run_subset_phylop.sh sub1
set -euo pipefail

TAG=${1:-sub1}
cd "$TAG"
mkdir -p mod2 log2 out2

phyloFit --tree tree/bang.treefile --subst-mod REV --nrates 4 \
  --precision MED --out-root mod2/concat aln/fit.fasta 2> log2/phylofit.err

phyloP --method LRT --mode CON --wig-scores --refidx 0 \
  mod2/concat.mod aln/concat.fasta \
  > out2/concat.CON.wig 2> log2/phylop.CON.err
