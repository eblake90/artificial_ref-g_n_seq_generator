# artificial_ref-g_n_seq_generator
This Python script is designed for bioinformatics simulations, specifically to generate mock genomic data. It can create a random genome sequence, introduce various types of genetic variants (such as duplications, insertions, deletions, and translocations), and simulate sequencing reads in FASTQ format.  Key features include:

**Random Genome Generation:** It can generate a genome sequence of specified length composed of random nucleotides (A, C, G, T).

**Variant Embedding:** The script is capable of randomly introducing different types of structural variants (like duplication, insertion, deletion, and translocation) into the generated genome.

**FASTQ Simulation:** It simulates sequencing reads by generating FASTQ entries from the modified genome. FASTQ format is a standard in bioinformatics for storing nucleotide sequences and their corresponding quality scores.

**Customizability and Reproducibility:** Users can specify the length of the genome, the number of variants, the length of sequencing reads, and the number of reads via command-line arguments. The script also supports setting a random seed for reproducibility, ensuring that the same inputs yield consistent outputs.

**Output Files:** The script writes the original genome, the variant-embedded genome, FASTQ entries for both normal and reverse complement sequences, and details of the introduced variants to separate files.

This script is useful for testing and developing bioinformatics tools and pipelines, as it allows for the creation of controlled, reproducible genomic datasets.


# Explaination of the Code
This Python script is designed to simulate the generation of genomic data for bioinformatics applications. The following are the main components and functionalities:

**Importing Libraries (Lines 1-3)**

argparse for parsing command-line arguments.

random for generating random numbers and sequences.

sys for system-specific parameters and functions (though it's imported but not used in this script).


## Function Definitions (Lines 5-63)

**generate_genome(length) (Lines 5-7):** Generates a random genome sequence of a specified length. It randomly chooses nucleotides ('A', 'C', 'G', 'T') to create a sequence.

**reverse_complement(seq) (Lines 9-12):** Creates the reverse complement of a DNA sequence. For each nucleotide, it finds the complement ('A' ↔ 'T', 'C' ↔ 'G') and reverses the order of the entire sequence.

**embed_variants(seq, num_variants) (Lines 14-48):** Introduces random genetic variants (like duplications, insertions, deletions, translocations) into the given sequence at random positions.

**generate_fastq(seq, read_length, num_reads) (Lines 50-58):** Simulates the generation of FASTQ entries. FASTQ format is used in bioinformatics to store sequence data along with quality scores. This function randomly selects subsequences (reads) from the given sequence and assigns random quality scores.

**write_to_file(filename, content) (Lines 60-63):** Writes the given content to a specified file.

## Main Function (Lines 65-85)

main() is the entry point of the script when executed. It uses the argparse library to parse command-line arguments for genome length (--genome_length), number of variants (--variants), read length (--read_length), and number of reads (--num_reads).

It first generates a random genome using generate_genome.

Then, it embeds variants into this genome using embed_variants.

Next, it generates FASTQ entries for both the normal and reverse complement of the variant-embedded genome using generate_fastq.

Finally, it writes the original genome, the FASTQ entries, and the variant details to separate files.
Script Execution Check (Line 87)

**if __name__ == '__main__': main():** This line checks if the script is being run as the main program and not being imported as a module in another script. If it's the main program, it calls the main() function.


