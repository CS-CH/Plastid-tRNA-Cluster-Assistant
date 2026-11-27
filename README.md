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

Single-nucleotide conservation of the four intron-containing plastid tRNAs (trnA-UGC, trnG-UCC, trnI-GAU and trnL-UAA) was estimated with phyloFit, phyloP and phastCons. Accessions were split into six disjoint subsets by make_species_subsets.py and run_subset_phylop.sh when a single all-accession PHAST run produced degenerate gamma-rate weights or matrix-exponential overflow during phyloFit.

After the subset analyses, `build_conservation_tracks.py` combines the scores on a common reference. `hierarchical_subset_tests.py` and `compute_subset_effects.py` compare introns with their flanking exons across subsets to estimate effect sizes and assess significance, rather than treating adjacent alignment sites as independent observations.

## 3. Intergenic-region analysis

`analyze_intergenic_regions.py` compares the intergenic regions of three tRNA gene pairs to identify sequence features associated with tRNA clustering. It summarizes length and base composition, identifies 8-mers enriched against a mononucleotide background, and evaluates their reliability using coverage, positional consistency, sequence complexity and possible flank-gene leakage.

High-confidence 8-mers are assembled into candidate core elements, while BLASTN is used to assess sequence similarity among regions. The resulting tables and figures can be used to examine conserved features that may contribute to tRNA clustering.

## Requirements

- Python 3 with NumPy, SciPy, openpyxl and matplotlib; Biopython is required by `add_sequences_to_gff.py`.
- tRNAscan-SE, PHAST (phyloFit, phyloP, phastCons), IQ-TREE, bedtools and NCBI BLAST+ on the command line.
