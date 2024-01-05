from Bio import SeqIO


def sliding_window(sequence, window_size=50, slide=1):
    counts = []
    for i in range(0, len(sequence) - window_size + 1, slide):
        window_sequence = sequence[i:i + window_size]
        g_count = window_sequence.count('G')
        counts.append(g_count)
    return counts


def average_g_content(sequence, window_size=50, slide=1):
    g_counts = sliding_window(sequence, window_size, slide)
    return sum(g_counts) / len(g_counts) if len(g_counts) > 0 else 0


def significant_g_content(average_g_contents, significance_level=0.05):
    sorted_averages = sorted(average_g_contents)
    index = int((1 - significance_level) * len(sorted_averages))
    cutoff = sorted_averages[index]
    return cutoff


def get_all_average_g_contents(input_file, num_sequences=1000):
    records = list(SeqIO.parse(input_file, "fasta"))
    shuffled_sequences = [str(record.seq) for record in records[:num_sequences]]

    average_g_contents = []
    for seq in shuffled_sequences:
        avg_g_content = average_g_content(seq)
        average_g_contents.append(avg_g_content)

    return average_g_contents


def calculate_threshold(input_file, num_sequences=1000):
    average_g_contents = get_all_average_g_contents(input_file, num_sequences)
    cutoff = significant_g_content(average_g_contents)
    return cutoff


def main():
    input_file = "shuffled_NC_045512.2.fasta"
    num_sequences = 1000
    average_g_contents = get_all_average_g_contents(input_file, num_sequences)

    cutoff = significant_g_content(average_g_contents)

    print(f"G-Content cutoff for P<0.05: {cutoff:.2f}")

    input("Press Enter to show average G contents.")

    print(average_g_contents)


if __name__ == "__main__":
    main()
