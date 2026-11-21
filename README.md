# Plastid tRNA Cluster Analysis Scripts

This repository provides a set of general-purpose scripts for analyzing plastid tRNAs that contain introns and gene clusters, together with the downstream conservation and intergenic-region analysis workflows.

## 1. Annotation and matrix generation

Prior to using these scripts, tRNA genes must be annotated directly from the genome using tRNAscan-SE. Similarity searches should be performed using BLASTN, Geneious, or other tools, followed by validation of the annotations using tRNAscan-SE. Verified hits should then be filtered and merged with the direct genome scan results to ensure accurate tRNA gene predictions.

After generating a combined annotation, the scripts `parse_trna_gff.py` and `add_sequences_to_gff.py` can be used to:

1. Calculate annotated exon regions and intergenic/spacer regions.
2. Extract corresponding sequences for downstream analyses.

Spacer and intergenic validation is critical before defining tRNA genes or tRNA gene clusters. Tools such as Geseq or Geneious can be used to examine both spacer regions between annotated exons of the same tRNA gene and intergenic regions between different tRNA genes. This step ensures that no other genes or genomic features are located within these regions. Manual curation is typically required to confirm introns and the true adjacency of specific tRNA genes.

Once exon connectivity and intergenic regions are confirmed, `tRNA_Matrix_Generator.py` can be used to generate three matrices for gene length, gene direction, and sequence information.

## 2. Conservation analysis of intron-containing tRNAs

`analysis/01_phylop_introns/` contains the per-nucleotide conservation analysis of the four intron-containing plastid tRNAs (trnA-UGC, trnG-UCC, trnI-GAU, trnL-UAA) across 3,699 plastome accessions. The accessions are split into six disjoint random subsets that are analysed independently and then combined, so that significance testing is based on independent estimates rather than on autocorrelated per-site scores.

Run in order:

1. `make_species_subsets.py` — split the accessions into K disjoint subsets, reusing the rows of the global alignment and pruning the full IQ-TREE topology (original branch lengths are kept and locked with `!` for phyloFit).
   Usage: `python3 make_species_subsets.py <subset 1..K> [K] [seed] [fit_cols]`
2. `run_subset_phylop.sh` — per subset: phyloFit (REV + gamma4) followed by phyloP (LRT, CON and ACC) and phastCons.
   Usage: `bash run_subset_phylop.sh sub1`
3. `build_conservation_tracks.py` — merge the six subset tracks on a common reference frame and write `perposition_<locus>.tsv` (per-nucleotide scores with 15/20/30 nt sliding-window means), `region_summary_allspecies.tsv`, `intron_occupancy.tsv` and `allspecies_track.npz`.
4. `hierarchical_subset_tests.py` — block-based (30 bp) and per-site intron-minus-flank contrasts per subset, followed by a one-sample t-test across the six subsets.
5. `compute_subset_effects.py` — per-subset effect sizes (intron minus flank, core columns with >=50% global coverage) and the corresponding cross-subset test.

## 3. Intergenic-region analysis

`analysis/02_intergenic_3pairs/analyze_intergenic_regions.py` analyses the intergenic regions of three tRNA pairs (trnG-UCC/trnR-UCU, trnI-GAU/trnA-UGC, trnL-UAA/trnF-GAA): length distributions, base composition, 8-mer scan with mononucleotide-model enrichment and reliability tiers, assembly of high-confidence 8-mers into core elements, and BLASTN all-against-all similarity. The G-fM pair is excluded (monocot-only subset without a tobacco record).

## Requirements

- Python 3 with numpy, openpyxl, matplotlib; scipy is optional (exact p-values); Biopython is required by `add_sequences_to_gff.py`.
- tRNAscan-SE, PHAST (phyloFit, phyloP, phastCons), IQ-TREE, bedtools and NCBI BLAST+ on the command line.
- The conservation workflow is designed for a Linux/WSL environment.

## Inputs and paths

The analysis scripts keep the working directories of the original analysis environment (`/root/trna/work`, `/mnt/c/Users/chang/Desktop/tRNA_migrate`, `C:\Users\chang\Desktop\tRNA_migrate`). Adjust these constants before rerunning.

Required input data:

- `chloALL/`: global alignment (`aln/concat.fasta`, `aln/layout.tsv`), `col_region.npz`, `regions.tsv`, `taxa.txt` and `tree/concat.treefile`.
- `records.tsv`: accession-to-group mapping for the intergenic analysis.
- `Dataset5_phylogeny.xlsx`: pair sheets containing the intergenic sequences.

See `analysis/README.md` for further details on the analysis workflows.
