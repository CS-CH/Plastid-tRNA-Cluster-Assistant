#!/usr/bin/env bash
set -euo pipefail

TAG=${1:-sub1}
cd "$TAG"
mkdir -p mod2 log2 out2 pc2

phyloFit --tree tree/bang.treefile --subst-mod REV --nrates 4 \
  --precision MED --out-root mod2/concat aln/fit.fasta 2> log2/phylofit.err

for mode in CON ACC; do
  phyloP --method LRT --mode "$mode" --wig-scores --refidx 0 \
    mod2/concat.mod aln/concat.fasta \
    > "out2/concat.$mode.wig" 2> "log2/phylop.$mode.err"
done

phastCons --rho 0.3 --target-coverage 0.10 --expected-length 30 \
  --most-conserved "pc2/el_${TAG}_main.bed" --score \
  aln/concat.fasta mod2/concat.mod \
  > "pc2/pp_${TAG}_main.wig" 2> log2/phastcons.err
