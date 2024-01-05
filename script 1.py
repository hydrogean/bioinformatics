from Bio import SeqIO
from Bio.Seq import Seq
import random

random.seed(42)


def shuffle_sequence(sequence):
    bases = list(sequence)
    random.shuffle(bases)
    return ''.join(bases)


input_file = "NC_045512.2.fasta"
output_file = "shuffled_NC_045512.2.fasta"
num_shuffle = 1000

original_record = SeqIO.read(input_file, "fasta")
original_sequence = str(original_record.seq)

shuffled_records = []
for i in range(num_shuffle):
    shuffled_sequence = shuffle_sequence(original_sequence)

    shuffled_record = SeqIO.SeqRecord(Seq(shuffled_sequence), id=f"shuffled_{i+1}", description="")

    shuffled_records.append(shuffled_record)

with open(output_file, "w") as output_handle:
    SeqIO.write(shuffled_records, output_handle, "fasta")

print(f"{num_shuffle} shuffled sequences written to {output_file}.")
