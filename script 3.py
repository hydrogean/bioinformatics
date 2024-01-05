from Bio import SeqIO
from script2 import calculate_threshold


def count_g_content(seq):
    return seq.count('G')


def sliding_window(seq, window_size=50, slide=1):
    for i in range(0, len(seq) - window_size + 1, slide):
        yield i, i + window_size, seq[i:i + window_size]


def main():
    fasta_file = "NC_045512.2.fasta"
    record = SeqIO.read(fasta_file, "fasta")
    sequence = str(record.seq)

    threshold = calculate_threshold("shuffled_NC_045512.2.fasta", num_sequences=1000)

    window_size = 50
    step = 1
    result = []

    for start, end, window_seq in sliding_window(sequence, window_size, step):
        g_content = count_g_content(window_seq)
        if g_content > threshold:
            result.append((record.id, start, end, g_content))

    output_file = "output_file.txt"
    with open(output_file, "w") as f:
        f.write("Name\tStart\tEnd\tG-Content\n")
        for name, start, end, g_content in result:
            f.write(f"{name}\t{start}\t{end}\t{g_content}\n")


if __name__ == "__main__":
    main()
