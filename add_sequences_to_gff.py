import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

path = ''

# Input/output files
gff_input = f'{path}/trna_annotations_final.gff3'
gff_output = f'{path}/trna_annotations_final_with_sequences.gff3'
fasta_input = f'{path}/Genome.fasta'
fasta_extracted = f'{path}/trna_extracted.fasta'

# Step 1: Extract sequences with bedtools
subprocess.run([
    'bedtools', 'getfasta',
    '-fi', fasta_input,
    '-bed', gff_input,
    '-fo', fasta_extracted
])

# Step 2: Load sequences
sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_extracted, 'fasta')}

# Step 3: Append sequences to GFF
output_lines = []
with open(gff_input, 'r') as gff:
    for line in gff:
        columns = line.strip().split('\t')
        start = int(columns[3]) - 1
        length = int(columns[4]) - start

        if length <= 15:
            seq_id = f'{columns[0]}:{start}-{columns[4]}'
            sequence = sequences.get(seq_id, '')

            if columns[6] == '+':
                columns.append(sequence)
            else:
                columns.append(str(Seq(sequence).reverse_complement()))

            output_lines.append('\t'.join(columns))

# Step 4: Write output
with open(gff_output, 'w') as out_file:
    out_file.write('\n'.join(output_lines))

print(f"[Done] Updated GFF with sequences written to: {gff_output}")
