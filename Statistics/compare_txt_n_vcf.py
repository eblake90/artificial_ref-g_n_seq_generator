# run in py3
# to run
# python compare_txt_n_vcf.py path_to_vcf_file.vcf path_to_txt_file.txt path_to_bam_file.bam path_to_fasta_file.fasta

import os
import datetime
import argparse
import pysam
import pandas as pd
from itertools import combinations



# Function to parse the VCF file and extract SV regions
# Function to parse the VCF file and extract SV regions
def parse_vcf(vcf_file):
    sv_regions = []
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            if 'SVTYPE' in record.info:
                sv_regions.append({
                    'chrom': record.chrom,  # Add chromosome info
                    'POS': record.pos,
                    'END': record.stop,
                    'SVTYPE': record.info['SVTYPE'],
                    'CIPOS': record.info.get('CIPOS', (0, 0)),
                    'CIEND': record.info.get('CIEND', (0, 0))
                })
    return sv_regions



# Function to extract sequence information using samtools
def extract_sequences(sv_regions, bam_file, fasta_file):
    sequences = []
    fasta = pysam.FastaFile(fasta_file)
    bam = pysam.AlignmentFile(bam_file)

    for sv in sv_regions:
        # Extracting the reference sequence
        ref_seq = fasta.fetch(reference=sv['chrom'], start=sv['POS'], end=sv['END'])
        # Extracting reads from the BAM file
        reads = bam.fetch(contig=sv['chrom'], start=sv['POS'], end=sv['END'])
        sequences.append({
            'SV': sv,
            'Sequence': ref_seq[:10] + '...',  # Store the sequence
            'Reads': list(reads)
        })
    return sequences



# Function to parse the TXT file for expected SVs
def parse_txt(txt_file):
    expected_svs = []
    with open(txt_file, 'r') as file:
        for line_number, line in enumerate(file, 1):
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            try:
                parts = line.split(',')
                sv_info = {p.split(':')[0].strip(): p.split(':')[1].strip() for p in parts}
                if 'Sequence' not in sv_info:  # Ensure 'Sequence' key exists
                    raise ValueError("Missing 'Sequence' in line")
                # Store only the first 10 base pairs of the sequence
                sv_info['Sequence'] = sv_info['Sequence'][:10] + '...'
                expected_svs.append(sv_info)
            except (IndexError, ValueError) as e:
                print(f"Error parsing line {line_number}: '{line}' - {e}")
    return expected_svs



# Function to calculate Hamming distance
def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


# Function to compare SVs and create the desired table
def compare_svs(predicted_svs, expected_svs):
    comparison_results = []

    for e_sv in expected_svs:
        e_sv_sequence = e_sv['Sequence']
        e_sv_type = e_sv['Type']  # Assumes 'Type' is a key in expected SVs
        similarities = []

        for p_sv in predicted_svs:
            p_sv_sequence = p_sv['Sequence']
            dist = hamming_distance(e_sv_sequence, p_sv_sequence)
            similarities.append((p_sv, dist))

        # Sort predicted SVs by similarity (distance)
        similarities.sort(key=lambda x: x[1])

        # Function to format top matches
        def format_match(match):
            return {
                'SV': match[0]['SV'],
                'Hamming_Distance': match[1],
                'First_10bp': match[0]['Sequence'][:10]
            }

        # Select top 3 most similar predicted SVs and format them
        top_matches = [format_match(sv) for sv in similarities[:3]]

        # Build the row for the comparison_results table
        row = {
            'Expected_SVs': e_sv,
            'Expected_SV_Type': e_sv_type
        }
        for i, match in enumerate(top_matches, start=1):
            row[f'Most_Probable_{i}'] = match.get('SV')
            row[f'Hamming_Distance_{i}'] = match.get('Hamming_Distance')
            row[f'First_10bp_{i}'] = match.get('First_10bp')

        comparison_results.append(row)

    return pd.DataFrame(comparison_results)

# Main function to parse arguments and run the script
def main():
    parser = argparse.ArgumentParser(
        description="Compare structural variants from a VCF file with expected SVs in a TXT file.")

    parser.add_argument("vcf_file", help="Path to the VCF file")
    parser.add_argument("txt_file", help="Path to the TXT file containing expected SVs")
    parser.add_argument("bam_file", help="Path to the BAM file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("--output", help="Path to save the output file. If not specified, the output will be printed.",
                        default="")

    args = parser.parse_args()

    sv_regions = parse_vcf(args.vcf_file)
    sequences = extract_sequences(sv_regions, args.bam_file, args.fasta_file)
    expected_svs = parse_txt(args.txt_file)
    comparison_table = compare_svs(sequences, expected_svs)


    if args.output:
        if not os.path.isdir(args.output):
            print(f"Error: The path '{args.output}' is not a directory.")
            return

        output_file = os.path.join(args.output, f"comparison_output.csv")

        comparison_table.to_csv(output_file, index=False)
        print(f"Output saved to {output_file}")
    else:
        print(comparison_table)

# Run the main function
if __name__ == "__main__":
    main()