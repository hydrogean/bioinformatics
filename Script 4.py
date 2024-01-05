import pandas as pd
from Bio import SeqIO


def calculate_g_content(sequence):
    g_count = sequence.upper().count("G")
    return g_count


mutations_data = pd.read_excel("mutations.xlsx")

fasta_file = "NC_045512.2.fasta"
genomic_record = SeqIO.read(fasta_file, "fasta")
genomic_sequence = genomic_record.seq

results = []

for index, row in mutations_data.iterrows():
    protein = row["Protein"]
    protein_change = row["ProteinChange"]
    genomic_location = row["GenomicLocation"]

    position = int(genomic_location)

    flanking_sequence = genomic_sequence[max(0, position-25):position+25]

    g_content = calculate_g_content(flanking_sequence)

    results.append({
        "Row": index + 2,
        "Protein": protein,
        "ProteinChange": protein_change,
        "Position": position,
        "GContent": g_content
    })

sorted_results = sorted(results, key=lambda x: x["GContent"], reverse=True)

output_file = "mutation_output.txt"

with open(output_file, 'w') as file:
    file.write("Top 10 Enriched Regions:\n")
    for region in sorted_results[:10]:
        file.write(f"Row: {region['Row']}, Gene: {region['Protein']}, Mutation: {region['ProteinChange']}, Position: {region['Position']}, G-content: {region['GContent']}\n")

print(f"Results written to {output_file}")
