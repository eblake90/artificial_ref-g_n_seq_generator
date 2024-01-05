# to run
# python fasta_n_fastq_generator.py --genome_length 1000000 --variants 100 --read_length 150 --num_reads 50000

import argparse
import random

def generate_genome(length):
    # Function to generate a random genome sequence of a specified length
    return ''.join(random.choice('ACGT') for _ in range(length))  # Generating random sequence of 'ACGT'

def reverse_complement(seq):
    # Function to generate the reverse complement of a DNA sequence
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # Mapping of nucleotides to their complements
    return ''.join(complement[base] for base in reversed(seq))  # Creating reverse complement sequence

def embed_variants(seq, num_variants):
    # Function to embed random variants into a sequence
    variants = []  # List to store variant details
    for _ in range(num_variants):
        variant_type = random.choice(['duplication', 'insertion', 'deletion', 'translocation'])  # Randomly selecting a variant type
        pos = random.randint(0, len(seq) - 1)  # Randomly selecting a position in the sequence

        # Handling duplication variants
        if variant_type == 'duplication':
            size = random.randint(1, 100)  # Determining the size of the duplication
            variant_seq = seq[pos:pos+size]  # Extracting the sequence to duplicate
            seq = seq[:pos] + variant_seq + seq[pos:]  # Inserting the duplicated sequence
            variants.append(f'Type: duplication, Location: {pos}, Length: {size}, Sequence: {variant_seq}')  # Recording the variant details

        # Handling insertion variants
        elif variant_type == 'insertion':
            size = random.randint(1, 100)  # Determining the size of the insertion
            insert_seq = ''.join(random.choice('ACGT') for _ in range(size))  # Generating a random sequence for insertion
            seq = seq[:pos] + insert_seq + seq[pos:]  # Inserting the new sequence
            variants.append(f'Type: insertion, Location: {pos}, Length: {size}, Sequence: {insert_seq}')  # Recording the variant details

        # Handling deletion variants
        elif variant_type == 'deletion':
            size = random.randint(1, 100)  # Determining the size of the deletion
            del_seq = seq[pos:pos+size]  # Extracting the sequence to be deleted
            seq = seq[:pos] + seq[pos+size:]  # Deleting the specified sequence
            variants.append(f'Type: deletion, Location: {pos}, Length: {size}, Sequence: {del_seq}')  # Recording the variant details

        # Handling translocation variants
        elif variant_type == 'translocation':
            size = random.randint(1, 100)  # Determining the size of the translocation
            end_pos = random.randint(0, len(seq) - size)  # Selecting a new position for the sequence
            trans_seq = seq[end_pos:end_pos+size]  # Extracting the sequence to translocate
            seq = seq[:pos] + trans_seq + seq[pos:]  # Inserting the translocated sequence
            seq = seq[:end_pos] + seq[end_pos+size:]  # Removing the original sequence
            variants.append(f'Type: translocation, From: {end_pos}, To: {pos}, Length: {size}, Sequence: {trans_seq}')  # Recording the variant details

    return seq, variants  # Returning the modified sequence and the list of variants

def generate_fastq(seq, read_length, num_reads):
    # Function to generate FASTQ entries
    fastq_entries = []  # List to store FASTQ entries
    for _ in range(num_reads):
        start = random.randint(0, len(seq) - read_length)  # Selecting a random start position for the read
        read_seq = seq[start:start + read_length]  # Extracting the sequence of the read
        quality_scores = ''.join([chr(random.randint(33, 73)) for _ in range(read_length)])  # Generating dummy quality scores
        fastq_entry = f'@read_{start}\n{read_seq}\n+\n{quality_scores}'  # Formatting the FASTQ entry
        fastq_entries.append(fastq_entry)  # Adding the entry to the list
    return fastq_entries  # Returning the list of FASTQ entries

def write_to_file(filename, content):
    # Function to write content to a file
    with open(filename, 'w') as file:
        file.write(content)  # Writing content to the specified file

def main():
    # Main function to handle workflow
    parser = argparse.ArgumentParser(description='Generate FASTA and FASTQ files with structural variants.')  # Creating an argument parser
    parser.add_argument('--genome_length', type=int, default=1000000)  # Adding argument for genome length
    parser.add_argument('--variants', type=int, default=100)  # Adding argument for number of variants
    parser.add_argument('--read_length', type=int, default=100)  # Adding argument for read length
    parser.add_argument('--num_reads', type=int, default=1000)  # Adding argument for number of reads
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')  # Adding argument for seed
    args = parser.parse_args()  # Parsing the arguments

    # Set the random seed for reproducibility
    random.seed(args.seed)  # This line sets the seed

    genome = generate_genome(args.genome_length)  # Generating a random genome

    genome_with_variants, variant_details = embed_variants(genome, args.variants)  # Embedding variants into the genome

    fastq_normal = generate_fastq(genome_with_variants, args.read_length, args.num_reads)  # Generating FASTQ for normal reads
    fastq_reverse_complement = generate_fastq(reverse_complement(genome_with_variants), args.read_length, args.num_reads)  # Generating FASTQ for reverse complement reads

    # Writing the generated data to files
    write_to_file('1e6-100v-150r-50000numr/reference_genome.fa', genome)
    write_to_file('1e6-100v-150r-50000numr/normal_reads.fastq', '\n'.join(fastq_normal))
    write_to_file('1e6-100v-150r-50000numr/reverse_complement_reads.fastq', '\n'.join(fastq_reverse_complement))
    write_to_file('1e6-100v-150r-50000numr/variant_details.txt', '\n'.join(variant_details))

if __name__ == '__main__':
    main()  # Running the main function if the script is executed directly