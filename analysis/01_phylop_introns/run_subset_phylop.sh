#!/bin/bash
# One species subset: phyloFit (fixed branch lengths) -> phyloP -> phastCons
# Usage: bash run_subset_phylop.sh sub1
set -u
TAG=${1:-sub1}
W=/root/trna/work/$TAG
cd "$W" || exit 1
mkdir -p mod2 log2 out2 pc2

S0=$(date +%s)
echo "[$(date +%T)] == $TAG start =="
phyloFit --tree tree/bang.treefile --subst-mod REV --nrates 4 \
         --precision MED --out-root mod2/concat aln/fit.fasta 2> log2/phylofit.err
RC=$?
echo "[$(date +%T)] phyloFit rc=$RC elapsed=$(( $(date +%s) - S0 ))s"
if [ $RC -ne 0 ]; then echo "FIT FAILED"; tail -3 log2/phylofit.err; exit 1; fi
grep -E "TRAINING_LNL|RATE_WEIGHTS|RATE_CONSTS" mod2/concat.mod | sed 's/^/  /'

for mode in CON ACC; do
  S1=$(date +%s)
  phyloP --method LRT --mode $mode --wig-scores --refidx 0 mod2/concat.mod aln/concat.fasta \
         > out2/concat.$mode.wig 2> log2/phylop.$mode.err
  RC=$?
  NV=$(grep -vc 'fixedStep\|track' out2/concat.$mode.wig)
  echo "[$(date +%T)] phyloP $mode rc=$RC values=$NV elapsed=$(( $(date +%s) - S1 ))s err=$(head -1 log2/phylop.$mode.err)"
done

S2=$(date +%s)
phastCons --rho 0.3 --target-coverage 0.10 --expected-length 30 \
          --most-conserved pc2/el_${TAG}_main.bed --score \
          aln/concat.fasta mod2/concat.mod > pc2/pp_${TAG}_main.wig 2> log2/phastcons.err
RC=$?
echo "[$(date +%T)] phastCons rc=$RC pp=$(wc -l < pc2/pp_${TAG}_main.wig) elements=$(wc -l < pc2/el_${TAG}_main.bed) elapsed=$(( $(date +%s) - S2 ))s"
echo "[$(date +%T)] == $TAG done, total $(( $(date +%s) - S0 ))s =="
