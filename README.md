# Plastid tRNA Cluster Analysis Scripts

This repository provides a set of general-purpose scripts for analyzing plastid tRNAs that contain introns and gene clusters.

Prior to using these scripts, tRNA genes must be annotated directly from the genome using tRNAscan-SE. Similarity searches should be performed using BLASTN, Geneious, or other tools, followed by validation of the annotations using tRNAscan-SE. Verified hits should then be filtered and merged with the direct genome scan results to ensure accurate tRNA gene predictions.

After generating a combined annotation, the scripts parse_trna_gff.py and add_sequences_to_gff.py can be used to:

Calculate annotated exon regions and intergenic/spacer regions.

Extract corresponding sequences for downstream analyses.

Spacer and intergenic validation is critical before defining tRNA genes or tRNA gene clusters. Tools such as Geseq or Geneious can be used to examine both spacer regions between annotated exons of the same tRNA gene and intergenic regions between different tRNA genes. This step ensures that no other genes or genomic features are located within these regions. Manual curation is typically required to confirm introns and the true adjacency of specific tRNA genes.

Once exon connectivity and intergenic regions are confirmed, tRNA_Matrix_Generator.py can be used to generate three matrices for gene length, gene direction, and sequence information.
