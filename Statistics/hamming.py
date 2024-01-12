# to run
# python script_name.py path_to_text_file_1.txt path_to_text_file_2.txt

import argparse
import os


import argparse
import os

def sliding_hamming_distance(short_seq, long_seq):
    """Calculate the Hamming distance using a sliding window approach."""
    min_distance = float('inf')
    for i in range(len(long_seq) - len(short_seq) + 1):
        window_seq = long_seq[i:i+len(short_seq)]
        distance = sum(ch1 != ch2 for ch1, ch2 in zip(short_seq, window_seq))
        if distance < min_distance:
            min_distance = distance
    return min_distance

def read_file_1(filepath):
    """Read data from file 1 and return sequences, ignoring NaN."""
    with open(filepath, 'r') as file:
        lines = file.readlines()
    return [line.strip().split(' // ') for line in lines[1:] if line.strip().split(' // ')[5] != 'NaN']

def read_file_2(filepath):
    """Read data from file 2 and return a list of sequences."""
    sequences = []
    with open(filepath, 'r') as file:
        for line in file:
            try:
                parts = line.strip().split(', ')
                sequences.append((parts[0].split(': ')[1], int(parts[2].split(': ')[1]), parts[3].split(': ')[1]))
            except IndexError:
                continue
    return sequences

def find_closest_sequences(file1_data, file2_sequences):
    """Find the closest sequence in file 2 for each sequence in file 1."""
    closest_sequences = []
    for data in file1_data:
        seq1 = data[5]
        min_distance = float('inf')
        closest_seq = None
        for seq2 in file2_sequences:
            if len(seq1) <= len(seq2[2]):
                distance = sliding_hamming_distance(seq1, seq2[2])
            else:
                distance = sliding_hamming_distance(seq2[2], seq1)
            if distance < min_distance:
                min_distance = distance
                closest_seq = seq2
        closest_sequences.append((data[0], data[4], closest_seq[0], closest_seq[1], min_distance))
    return closest_sequences

def write_output(filepath, data):
    """Write the output to a file."""
    with open(filepath, 'w') as file:
        file.write("ID // SVTYPE // Type / Length // Hamming Distance\n")
        for item in data:
            file.write(f"{item[0]} // {item[1]} // {item[2]} / {item[3]} // {item[4]}\n")

def main(file1_path, file2_path):
    file1_data = read_file_1(file1_path)
    file2_sequences = read_file_2(file2_path)
    closest_sequences = find_closest_sequences(file1_data, file2_sequences)

    output_filename = os.path.splitext(os.path.basename(file2_path))[0] + "_hammed_svs.txt"
    output_path = os.path.join(os.path.dirname(file2_path), output_filename)
    write_output(output_path, closest_sequences)
    print(f"Output saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find sequences from file 1 with the lowest Hamming distance to each sequence in file 2 using a sliding window approach.")
    parser.add_argument("file1_path", help="Path to the first text file")
    parser.add_argument("file2_path", help="Path to the second text file")
    args = parser.parse_args()

    main(args.file1_path, args.file2_path)
