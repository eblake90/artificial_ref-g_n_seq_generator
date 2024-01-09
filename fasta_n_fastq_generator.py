# run in py3

import argparse
import random
import os

def generate_genome(length):
    # Function to generate a random genome sequence of a specified length
    return ''.join(random.choice('ACGT') for _ in range(length))  # Generating random sequence of 'ACGT'

def reverse_complement(seq):
    # Function to generate the reverse complement of a DNA sequence
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # Mapping of nucleotides to their complements
    return ''.join(complement.get(base, base) for base in reversed(seq))  # Creating reverse complement sequence

def embed_variants(seq, num_variants):
    # Function to embed random variants into a sequence
    variants = []  # List to store variant details
    for _ in range(num_variants):
        variant_type = random.choice(['duplication', 'insertion', 'deletion', 'translocation'])  # Randomly selecting a variant type
        pos = random.randint(0, len(seq) - 1)  # Randomly selecting a position in the sequence

        # Handling duplication variants
        if variant_type == 'duplication':
            size = random.randint(50, 100)  # Determining the size of the duplication
            variant_seq = seq[pos:pos+size]  # Extracting the sequence to duplicate
            seq = seq[:pos] + variant_seq + seq[pos:]  # Inserting the duplicated sequence
            variants.append(f'Type: duplication, Location: {pos}, Length: {size}, Sequence: {variant_seq}')  # Recording the variant details

        # Handling insertion variants
        elif variant_type == 'insertion':
            size = random.randint(50, 100)  # Determining the size of the insertion
            insert_seq = ''.join(random.choice('ACGT') for _ in range(size))  # Generating a random sequence for insertion
            seq = seq[:pos] + insert_seq + seq[pos:]  # Inserting the new sequence
            variants.append(f'Type: insertion, Location: {pos}, Length: {size}, Sequence: {insert_seq}')  # Recording the variant details

        # Handling deletion variants
        elif variant_type == 'deletion':
            size = random.randint(50, 100)  # Determining the size of the deletion
            del_seq = seq[pos:pos+size]  # Extracting the sequence to be deleted
            seq = seq[:pos] + seq[pos+size:]  # Deleting the specified sequence
            variants.append(f'Type: deletion, Location: {pos}, Length: {size}, Sequence: {del_seq}')  # Recording the variant details

        # Handling translocation variants
        elif variant_type == 'translocation':
            size = random.randint(50, 100)  # Determining the size of the translocation
            end_pos = random.randint(0, len(seq) - size)  # Selecting a new position for the sequence
            trans_seq = seq[end_pos:end_pos+size]  # Extracting the sequence to translocate
            seq = seq[:pos] + trans_seq + seq[pos:]  # Inserting the translocated sequence
            seq = seq[:end_pos] + seq[end_pos+size:]  # Removing the original sequence
            variants.append(f'Type: translocation, From: {end_pos}, To: {pos}, Length: {size}, Sequence: {trans_seq}')  # Recording the variant details

    return seq, variants  # Returning the modified sequence and the list of variants

def generate_fragments(seq, fragment_length, num_fragments):
    fragments = []# List to store FASTQ entries
    for _ in range(num_fragments):
        start = random.randint(0, len(seq) - fragment_length) # Selecting a random start position for the read
        fragment = seq[start:start + fragment_length]  # Extracting the sequence of the read
        fragments.append(fragment) # Adding the entry to the list
    return fragments  # Returning the list of FASTQ entries

def generate_paired_end_fastq(fragments, read_length):
    paired_fastq_entries = []
    for fragment in fragments:
        if len(fragment) < 2 * read_length:
            continue

        forward_read = fragment[:read_length]
        forward_quality_scores = ''.join([chr(random.randint(33, 73)) for _ in range(read_length)])

        reverse_read = reverse_complement(fragment[-read_length:])
        reverse_quality_scores = ''.join([chr(random.randint(33, 73)) for _ in range(read_length)])

        paired_fastq_entries.append((forward_read, forward_quality_scores, reverse_read, reverse_quality_scores))
    return paired_fastq_entries

def write_reads_to_separate_files(directory, forward_filename, reverse_filename, paired_fastq_entries):
    with open(os.path.join(directory, forward_filename), 'w') as forward_file, \
         open(os.path.join(directory, reverse_filename), 'w') as reverse_file:
        for i, entry in enumerate(paired_fastq_entries):
            read_name_base = f"@read_{i}"  # Common base name for the paired read
            forward_file.write(f"{read_name_base}/1\n{entry[0]}\n+\n{entry[1]}\n")
            reverse_file.write(f"{read_name_base}/2\n{entry[2]}\n+\n{entry[3]}\n")

def write_to_file(directory, filename, header, content):
    with open(os.path.join(directory, filename), 'w') as file:
        file.write(header + '\n' + content)

def main():
    # Main function to handle workflow
    parser = argparse.ArgumentParser(description='Generate FASTA and FASTQ files with structural variants.')  # Creating an argument parser
    parser.add_argument('--genome_using', choices=['real', 'fake'], required=True, help='Use real or artificially generated genome (fake)')
    parser.add_argument('--real_genome_path', type=str,
                        help='Path to the real genome file (required if --genome_using real)')
    parser.add_argument('--genome_length', type=int, default=1000000)  # Adding argument for genome length
    parser.add_argument('--variants', type=int, default=100)  # Adding argument for number of variants
    parser.add_argument('--read_length', type=int, default=100)  # Adding argument for read length
    parser.add_argument('--num_reads', type=int, default=1000)  # Adding argument for number of reads
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')  # Adding argument for seed
    parser.add_argument('--output_dir', required=True, help='Output directory for saving files')

    args = parser.parse_args()  # Parsing the arguments

    # Validate inputs
    if args.genome_using == 'real' and not args.real_genome_path:
        parser.error("--real_genome_path is required when --genome_using is set to 'real'")

    # Set the random seed for reproducibility
    random.seed(args.seed)  # This line sets the seed

    if args.genome_using == 'fake':
        # Generate artificial genome
        genome = generate_genome(args.genome_length) # Generating a random genome
    else:
        # Use real genome
        with open(args.real_genome_path, 'r') as file:
            genome = file.read()

    # Generating the genome with variants
    genome_with_variants, variant_details = embed_variants(genome, args.variants)  # Embedding variants into the genome

    # Preparing for paired-end read generation
    fragment_length = 3 * args.read_length  # Adjust as needed
    num_fragments = args.num_reads // 2  # Assuming two reads per fragment #########

    fragments = generate_fragments(genome_with_variants, fragment_length, num_fragments)
    paired_end_reads = generate_paired_end_fastq(fragments, args.read_length)

    # Constructing file names based on input arguments
    if args.genome_using == 'fake':
        base_filename = f"output_fake_{args.genome_length}_{args.variants}v_{args.read_length}rl_{args.num_reads}rn_{args.seed}s"
    else:
        base_filename = f"output_real_{args.variants}v_{args.read_length}rl_{args.num_reads}_{args.seed}s"

    output_subfolder = os.path.join(args.output_dir, f"{base_filename}")
    os.makedirs(output_subfolder, exist_ok=True)

    # Constructing the header for the reference genome
    genome_header = f">ref_genome_{base_filename}"

    # Writing the generated data to files in the specified directory
    if args.genome_using == 'fake':
        write_to_file(output_subfolder, f'reference_genome_{base_filename}.fa', genome_header, genome)
    write_to_file(output_subfolder, f'variant_details_{base_filename}.txt', '', '\n'.join(variant_details))
    write_reads_to_separate_files(output_subfolder, f'forward_reads_{base_filename}.fastq', f'reverse_reads_{base_filename}.fastq', paired_end_reads)

if __name__ == '__main__':
    main()  # Running the main function if the script is executed directly