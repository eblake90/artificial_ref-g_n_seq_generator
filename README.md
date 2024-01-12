# artificial_ref-g_n_seq_generator Description
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

---

## `fasta_n_fastq_generator.py`
Usage Instructions
This script can be used in two modes: using a fake (artificially generated) genome or using a real (user-provided) genome.

#### Using an Artificial Genome
To run the script with an artificially generated genome, use the following command:
``` shell 
python fasta_n_fastq_generator.py --genome_using fake --genome_length 1000000 --variants 100 --read_length 150 --num_reads 50000 --seed 42 --output_dir /path/to/output
```
##### Parameters:
**--genome_using fake:** Indicates the use of an artificial genome.

**--genome_length:** The length of the artificial genome to be generated.

Other parameters (--variants, --read_length, --num_reads, --seed, --output_dir) are used for variant embedding, read generation, and output file management.

#### Using a Real Genome
To run the script with a real genome provided by the user, use the following command:
```shell
python fasta_n_fastq_generator.py --genome_using real --real_genome_path /my/path/ref.fa --variants 100 --read_length 150 --num_reads 50000 --seed 42 --output_dir /path/to/output
```
##### Parameters:

**--genome_using real:** Indicates the use of a real genome.

**--real_genome_path:** Path to the real genome file.

Other parameters (--variants, --read_length, --num_reads, --seed, --output_dir) are used similarly as in the fake genome mode, but the real genome is used as the base sequence.

Note: In the 'real' mode, the --genome_length parameter is not required as the length will be determined by the provided genome file.

---

## `fasta_to_fa.py`

#### Description
This script is used to convert a FASTA file to a FA file. It reads a given FASTA file and writes its content into an FA file without making any modifications to the content.

#### Usage
```bash
python fasta_to_fa.py <input.fasta> <output.fa>
```

#### Parameters
- `input.fasta`: The input file in FASTA format.
- `output.fa`: The desired output file in FA format.

#### Example
```bash
python fasta_to_fa.py input.fasta output.fa
```

---

## `fastq_to_bam.py`

#### Description
This script processes sequencing reads by aligning them to a reference genome using BWA, converting the resulting SAM file to a BAM file using Samtools, sorting the BAM file, and indexing it for further analysis.

#### Usage
```bash
python fastq_to_bam.py --ref_genome <my.fa> --fastq1 <my1.fastq> --fastq2 <my2.fastq> --output </myoutput>
```

#### Parameters
- `--ref_genome`: The reference genome file in FASTA format.
- `--fastq1`: The first FASTQ file containing sequencing reads.
- `--fastq2`: The second FASTQ file containing sequencing reads.
- `--output`: Prefix for the output files (e.g., `/myoutput` for `/myoutput.sam`, `/myoutput.bam`, etc.).

#### Example
```bash
python fastq_to_bam.py --ref_genome my.fa --fastq1 my1.fastq --fastq2 my2.fastq --output /myoutput
```

---

## `indexing_ref_genome.py`

#### Description
This script is designed for indexing a reference genome. It creates a BWA index of the reference genome and a .fa.fai index file using Samtools. These indices are essential for various bioinformatics analyses that involve the reference genome.

#### Usage
```bash
python indexing_ref_genome.py --ref_genome <ref_genome.fa>
```

#### Parameters
- `--ref_genome`: The reference genome file in FASTA format that needs to be indexed.

#### Example
```bash
python indexing_ref_genome.py --ref_genome ref_genome.fa
```

---

# SAM File Read Summary

This section summarizes key information from the SAM file generated from aligning FASTQ files.

## Table: Alignment Summary

| Read Name | Flag | Reference Sequence Name                                      | Position | Mapping Quality | CIGAR String | RNEXT | PNEXT  | TLEN | Alignment Score (AS) | Number of Mismatches (NM) |
|-----------|------|--------------------------------------------------------------|----------|-----------------|--------------|-------|--------|------|----------------------|---------------------------|
| read_0    | 99   | ref_genome_output_fake_1000000_25v_150rl_50000rn_42s         | 475795   | 60              | 150M         | =     | 476095 | 450  | 150                  | 0                         |
| read_0    | 147  | ref_genome_output_fake_1000000_25v_150rl_50000rn_42s         | 476095   | 60              | 150M         | =     | 475795 | -450 | 150                  | 0                         |

## Explanation of Columns

- **Read Name**: Identifier of the read.
- **Flag**: Bitwise flag (99 and 147 indicate proper pairing for paired-end reads).
- **Reference Sequence Name**: The name of the reference sequence where the read aligns.
- **Position**: The starting position of the alignment on the reference sequence.
- **Mapping Quality**: The quality score of the alignment (60 is a high score, indicating high confidence).
- **CIGAR String**: Encodes the alignment (150M means 150 bases matched).
- **RNEXT**: Reference name of the mate/next read (‘=’ indicates the same chromosome).
- **PNEXT**: Position of the mate/next read.
- **TLEN**: Observed Template Length (450 and -450 indicate the size of the DNA fragment; negative value for the reverse strand).
- **Alignment Score (AS)**: A score indicating how well the read aligns to the reference.
- **Number of Mismatches (NM)**: The number of mismatches in the alignment.

### Significance of Flags 99 and 147
In the SAM file, the flags 99 and 147 have specific meanings related to paired-end sequencing:

**Flag 99**: Indicates that the read is part of a read pair (paired-end), it's the first read in the pair (read_1), the read is mapped in a proper pair, and it's mapped to the forward strand.

**Flag 147**: Indicates that the read is part of a read pair, it's the second read in the pair (read_2), the read is mapped in a proper pair, and it's mapped to the reverse strand.

These flags are a bitwise representation of various properties of each read. For example:

**99** in binary is **1100011**, indicating the read is paired, mapped in a proper pair, not unmapped, the first in pair, and not on the reverse strand.

**147** in binary is **10010011**, indicating the read is paired, mapped in a proper pair, not unmapped, the second in pair, and on the reverse strand.

---

