import pandas as pd
from pyfaidx import Fasta
import argparse

def extract_sequence(fasta, chrom, start, end):
    if chrom not in fasta:
        raise ValueError(f"Chromosome {chrom} not found in the reference genome.")

    full_seq = str(fasta[chrom][start-1:end])
    max_length = 5000  # Maximum sequence length to extract
    return full_seq[:max_length]

def extract_svtype(info):
    for item in info.split(';'):
        if item.startswith("SVTYPE"):
            return item.split('=')[1]
    return "Unknown"

def process_csv(csv_file, fa_file, output_file):
    df = pd.read_csv(csv_file)
    fasta = Fasta(fa_file)
    output = ["ID // POS // END // SVLEN // SVTYPE // SEQUENCE"]

    for _, row in df.iterrows():
        svtype = extract_svtype(row['INFO'])
        if svtype == "DEL":
            seq = extract_sequence(fasta, row['CHROM'], row['POS'], int(row['INFO'].split(';')[0].split('=')[1]))
        else:
            seq = "NaN"

        output.append(f"{row['ID']} // {row['POS']} // {row['INFO'].split(';')[0].split('=')[1]} // {row['INFO'].split(';')[2].split('=')[1]} // {svtype} // {seq}")

    with open(output_file, "w") as f:
        for line in output:
            f.write(line + "\n")

def main():
    parser = argparse.ArgumentParser(description='Process a CSV file with genome data and extract sequences.')
    parser.add_argument('csv_file', help='CSV file containing the genome data')
    parser.add_argument('fa_file', help='Reference genome file (FASTA)')
    parser.add_argument('output', help='Output file path for the extracted sequences')
    args = parser.parse_args()
    process_csv(args.csv_file, args.fa_file, args.output)

if __name__ == "__main__":
    main()
