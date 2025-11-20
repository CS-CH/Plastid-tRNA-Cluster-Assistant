import re
from typing import Dict, List, Any, Optional

# ------------------------------
# File paths and parameters
# ------------------------------
path = ""  # Fill with directory if needed, e.g., "data/"
chrom_size_file = f"{path}chrom_sizes.tsv"          # Chromosome sizes file: tab-separated, chrom_name <tab> size
gff_raw_file = f"{path}trna_annotations_raw.gff3"  # Raw tRNA annotation GFF3 file
gff_final_file = f"{path}trna_annotations_final.gff3"  # Output GFF file with gene and distance entries

# ------------------------------
# Gene order explanation
# ------------------------------
# gene_order string defines the order of tRNA genes/features for distance calculation.
# Underscore "_" separates adjacent genes in order.
# Slash "/" separates different exons/features of the same gene.
# The number prefix in gene name (e.g., 2trnG-UCC) indicates the total number of exons:
#   - "2trnG-UCC" means tRNA trnG-UCC has 2 exons (contains intron)
#   - "1trnX-XXX" means tRNA has only 1 exon (no intron)
# Example:
# "2trnG-UCC/1_2trnG-UCC/2_1trnR-UCU_2trnL-UAA/1_2trnL-UAA/2"
# Interpretation:
# - "2trnG-UCC/1" : first exon of tRNA trnG-UCC (which has 2 exons total)
# - "2trnG-UCC/2" : second exon of tRNA trnG-UCC
# - "1trnR-UCU"   : single-exon tRNA trnR-UCU (no intron)
# "_" indicates sequential order for calculating minimum distances between adjacent genes/features

gene_order = (
    "2trnG-UCC/1_2trnG-UCC/2_1trnR-UCU_2trnL-UAA/1_2trnL-UAA/2_"
    "1trnF-GAA_2trnI-GAU/1_2trnI-GAU/2_2trnA-UGC/1_2trnA-UGC/2"
)

# ------------------------------
# Function: parse_chromosome_sizes
# ------------------------------
def parse_chromosome_sizes(genome_file: str) -> Dict[str, int]:
    """
    Read chromosome size file and return chromosome lengths.

    Parameters
    ----------
    genome_file : str
        Path to chromosome size file, tab-separated: chrom_name <tab> size

    Returns
    -------
    Dict[str, int]
        Mapping of chromosome name to its length
    """
    chromosome_sizes: Dict[str, int] = {}
    with open(genome_file, "r") as f:
        for line in f:
            chrom_name, chrom_size = line.strip().split("\t")
            chromosome_sizes[chrom_name] = int(chrom_size)
    return chromosome_sizes

# ------------------------------
# Function: parse_gff
# ------------------------------
def parse_gff(gff_file: str) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    """
    Parse tRNA GFF file and store gene information by chromosome and gene name.

    Returns structure:
        {
            chromosome_name: {
                gene_name: [
                    {"start": int, "end": int, "strand": "+/-", "feature_number": str},
                ],
            },
        }
    """
    chrom_gene_positions: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}

    with open(gff_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue  # Skip malformed lines

            attributes_raw = cols[8].strip()
            gene_name_match = re.search(r"Name=([^;]+)", attributes_raw)
            if not gene_name_match:
                continue

            gene_name = gene_name_match.group(1)  # Gene name, e.g., "2trnG-UCC"
            chrom_name = cols[0]                  # Chromosome name
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]                      # "+" or "-"
            feature_number = attributes_raw.split("/")[-2]  # Exon or feature number, e.g., "1" or "2"

            gene_info = {
                "start": start,
                "end": end,
                "strand": strand,
                "feature_number": feature_number,
            }

            # Initialize dictionary and append gene info
            chrom_gene_positions.setdefault(chrom_name, {}).setdefault(gene_name, []).append(gene_info)

    return chrom_gene_positions

# ------------------------------
# Function: calculate_min_gene_distances
# ------------------------------
def calculate_min_gene_distances(
    gene_order: str,
    chrom_gene_positions: Dict[str, Dict[str, List[Dict[str, Any]]]],
    chromosome_sizes: Dict[str, int]
) -> List[Dict[str, Any]]:
    """
    Calculate minimum distances between adjacent genes/features according to a given gene_order.

    Strand-sensitive and wrap-around logic is retained:
        - Same strand: calculate distance normally; wrap-around considered for circular chromosomes
        - Opposite strand: assign a large value (ignored)

    Returns
    -------
    List[Dict[str, Any]]
        Each dictionary contains:
        {
            'chromosome_name': str,
            'geneA': str or None,
            'geneB': str or None,
            'distance': int or None,
            'feature_numberA': str or None,
            'feature_numberB': str or None,
            'positionA': dict or None,
            'positionB': dict or None
        }
    """
    ordered_genes = gene_order.split('_')
    num_genes = len(ordered_genes)
    min_distances: List[Dict[str, Any]] = []

    def calculate_distance(
        geneA: str,
        geneB: str,
        chromosome_size: int,
        chromosome_name: str
    ) -> tuple[int, Optional[str], Optional[str], Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:

        min_distance = float('inf')
        min_featureA: Optional[str] = None
        min_featureB: Optional[str] = None
        min_positionA: Optional[Dict[str, Any]] = None
        min_positionB: Optional[Dict[str, Any]] = None

        gene_positions = chrom_gene_positions[chromosome_name]

        for posA in gene_positions[geneA]:
            startA, endA, strandA = posA['start'], posA['end'], posA['strand']
            for posB in gene_positions[geneB]:
                startB, endB, strandB = posB['start'], posB['end'], posB['strand']

                # ---------- Original logic: strand-sensitive + wrap-around ----------
                if strandA == strandB:
                    if strandA == '+':
                        if startB > endA:
                            current_distance = startB - endA - 1
                        else:
                            current_distance = chromosome_size - endA + startB - 1
                    else:  # strandA == '-'
                        if endB < startA:
                            current_distance = startA - endB - 1
                        else:
                            current_distance = chromosome_size - endB + startA - 1
                else:
                    current_distance = 1e10  # Ignore opposite strand

                if current_distance < min_distance:
                    min_distance = current_distance
                    min_featureA = posA['feature_number']
                    min_featureB = posB['feature_number']
                    min_positionA = posA
                    min_positionB = posB

        return min_distance, min_featureA, min_featureB, min_positionA, min_positionB

    # Iterate chromosomes and gene order
    for chrom_name, gene_dict in chrom_gene_positions.items():
        for i in range(num_genes - 1):
            geneA, geneB = ordered_genes[i], ordered_genes[i + 1]

            if geneA in gene_dict and geneB in gene_dict:
                dist, featureA, featureB, posA, posB = calculate_distance(
                    geneA, geneB, chromosome_sizes[chrom_name], chrom_name
                )
                min_distances.append({
                    'chromosome_name': chrom_name,
                    'geneA': geneA,
                    'geneB': geneB,
                    'distance': dist,
                    'feature_numberA': featureA,
                    'feature_numberB': featureB,
                    'positionA': posA,
                    'positionB': posB
                })
            elif geneA in gene_dict:
                posA = gene_dict[geneA][0]
                min_distances.append({
                    'chromosome_name': chrom_name,
                    'geneA': geneA,
                    'geneB': None,
                    'distance': None,
                    'feature_numberA': posA['feature_number'],
                    'feature_numberB': None,
                    'positionA': posA,
                    'positionB': None
                })
            elif geneB in gene_dict:
                posB = gene_dict[geneB][0]
                min_distances.append({
                    'chromosome_name': chrom_name,
                    'geneA': None,
                    'geneB': geneB,
                    'distance': None,
                    'feature_numberA': None,
                    'feature_numberB': posB['feature_number'],
                    'positionA': None,
                    'positionB': posB
                })

    return min_distances

# ------------------------------
# Function: generate_output_gff
# ------------------------------
def generate_output_gff(output_filename: str, min_distances: List[Dict[str, Any]]) -> None:
    """
    Generate GFF file including tRNA genes and distances between gene pairs.
    """
    with open(output_filename, 'w') as f:
        for info in min_distances:
            chrom = info['chromosome_name']
            geneA, geneB = info.get('geneA'), info.get('geneB')
            distance = info.get('distance')
            posA, posB = info.get('positionA'), info.get('positionB')
            featureA, featureB = info.get('feature_numberA'), info.get('feature_numberB')

            if posA:
                f.write(f"{chrom}\t.\ttRNA\t{posA['start']}\t{posA['end']}\t.\t{posA['strand']}\t.\tName={geneA};Feature={featureA}\n")
            if posB:
                f.write(f"{chrom}\t.\ttRNA\t{posB['start']}\t{posB['end']}\t.\t{posB['strand']}\t.\tName={geneB};Feature={featureB}\n")

            if posA and posB and distance is not None:
                if posB['strand'] == '+':
                    start_dist, end_dist = posA['end'] + 1, posB['start'] - 1
                else:
                    start_dist, end_dist = posB['end'] + 1, posA['start'] - 1
                f.write(f"{chrom}\t.\tdistance\t{start_dist}\t{end_dist}\t.\t{posA['strand']}_{posB['strand']}\t.\tGenePair={geneA}_{geneB};Distance={distance}\n")

# ------------------------------
# Function: process_output_gff
# ------------------------------
def process_output_gff(output_filename: str) -> None:
    """
    Post-process the GFF file:
    - Remove duplicate lines
    - Correct start/end coordinates
    - Mark features with length 0 or 1 as Null
    """
    lines_seen = set()
    processed_lines = []

    with open(output_filename, 'r') as f:
        for line in f:
            if line not in lines_seen:
                lines_seen.add(line)
                processed_lines.append(line)

    with open(output_filename, 'w') as f:
        for line in processed_lines:
            cols = line.strip().split('\t')
            if len(cols) >= 5:
                start, end = int(cols[3]), int(cols[4])
                if abs(end - start) <= 1:
                    cols[0] = f"{cols[0]}/Null"
                elif start > end:
                    cols[3], cols[4] = cols[4], cols[3]
            f.write('\t'.join(cols) + '\n')

# ------------------------------
# Main workflow
# ------------------------------
chromosome_sizes = parse_chromosome_sizes(chrom_size_file)
gene_positions = parse_gff(gff_raw_file)
min_distances = calculate_min_gene_distances(gene_order, gene_positions, chromosome_sizes)
generate_output_gff(gff_final_file, min_distances)
process_output_gff(gff_final_file)

print(f"GFF file processing completed. Result saved to: {gff_final_file}")
