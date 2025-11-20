path = ''

# ------------------------------
# File paths
# ------------------------------
gff_file_path = f'{path}/trna_annotations_final_with_sequences.gff3'
length_output_path = f'{path}/trna_gene_length.matrix'
direction_output_path = f'{path}/trna_gene_direction.matrix'
sequence_output_path = f'{path}/trna_gene_sequence.matrix'

# ------------------------------
# Gene order string
# ------------------------------
gene_order = '2trnG-UCC/1_2trnG-UCC/2_1trnR-UCU_2trnL-UAA/1_2trnL-UAA/2_1trnF-GAA_2trnI-GAU/1_2trnI-GAU/2_2trnA-UGC/1_2trnA-UGC/2'

# ------------------------------
# Function to calculate gene length
# ------------------------------
def calculate_gene_length(start, end):
    return int(end) - int(start) + 1

# ------------------------------
# Function to create gene list including adjacent pairs
# ------------------------------
def gene_list(gene_order):
    gene_pair_list = gene_order.split('_')
    gene_order_list = [gene_pair_list[i] + '_' + gene_pair_list[i+1] for i in range(0, len(gene_pair_list)-1, 1)]

    merged_list = []
    max_length = max(len(gene_order_list), len(gene_pair_list))
    for i in range(max_length):
        if i < len(gene_pair_list):
            merged_list.append(gene_pair_list[i])
        if i < len(gene_order_list):
            merged_list.append(gene_order_list[i])
    return merged_list

# ------------------------------
# Function to process GFF file
# ------------------------------
def process_gff(file_path, genes):
    result_len = {}
    result_dir = {}
    result_seq = {}

    with open(file_path, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if fields[0].split('/')[-1] != 'Null':
                    chrom, direction, attributes, sequence = fields[0], fields[6], fields[8], fields[9]
                    parsed_attributes = fields[8].split('=')[1].split(';')[0]

                    for gene in genes:
                        if gene == parsed_attributes:
                            if '_' not in gene:
                                gene_length = calculate_gene_length(fields[3], fields[4])
                            else:
                                gene_length = fields[8].split('=')[-1]
                                if gene_length == '1e+1':
                                    gene_length = '--'
                            if chrom not in result_len:
                                result_len[chrom] = {}
                                result_dir[chrom] = {}
                                result_seq[chrom] = {}
                            result_len[chrom][gene] = gene_length
                            result_dir[chrom][gene] = direction
                            result_seq[chrom][gene] = sequence

    return result_len, result_dir, result_seq

# ------------------------------
# Main workflow
# ------------------------------
def main(gene_order, gff_file_path, length_output_path, direction_output_path, sequence_output_path):
    genes = gene_list(gene_order)
    gene_lengths, gene_direction, gene_sequences = process_gff(gff_file_path, genes)

    # Output gene lengths
    with open(length_output_path, 'w') as len_file:
        len_file.write('\t'.join(['Chromosome'] + genes) + '\n')
        for chrom, gene_length_dict in gene_lengths.items():
            lengths = [str(gene_length_dict[gene]) if gene in gene_length_dict else '' for gene in genes]
            len_file.write('\t'.join([chrom] + lengths) + '\n')

    # Output gene directions
    with open(direction_output_path, 'w') as dir_file:
        dir_file.write('\t'.join(['Chromosome'] + genes) + '\n')
        for chrom, gene_direction_dict in gene_direction.items():
            direction = [str(gene_direction_dict[gene]) if gene in gene_direction_dict else '' for gene in genes]
            dir_file.write('\t'.join([chrom] + direction) + '\n')

    # Output gene sequences
    with open(sequence_output_path, 'w') as seq_file:
        seq_file.write('\t'.join(['Chromosome'] + genes) + '\n')
        for chrom, gene_sequence_dict in gene_sequences.items():
            sequence = [str(gene_sequence_dict[gene]) if gene in gene_sequence_dict else '' for gene in genes]
            seq_file.write('\t'.join([chrom] + sequence) + '\n')


if __name__ == '__main__':
    main(gene_order, gff_file_path, length_output_path, direction_output_path, sequence_output_path)

print(f'GFF file has been processed. Results saved as {length_output_path}, {direction_output_path}, {sequence_output_path}')
