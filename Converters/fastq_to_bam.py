# run on py3
# to run
# python fastq_to_bam.py --ref_genome my.fa --fastq1 my1.fastq --fastq2 my2.fastq --output /myoutput
#python /home/eblake/Documents/manta/SV_finder_manta/results/fastq_to_bam.py --ref_genome /home/eblake/Documents/making_dna/output_360000000_100v_150rl_50000rn_42s/reference_genome.fa --fastq1 /home/eblake/Documents/making_dna/output_360000000_100v_150rl_50000rn_42s/forward_reads.fastq --fastq2 /home/eblake/Documents/making_dna/output_360000000_100v_150rl_50000rn_42s/reverse_reads.fastq --output /home/eblake/Documents/github/1.1/artificial_ref-g_n_seq_generator/artificial_stuff

import subprocess
import argparse

def run_command(command, output_file=None):
    """ Run a shell command """
    try:
        if output_file:
            with open(output_file, 'w') as f:
                subprocess.run(command, check=True, stdout=f)
        else:
            subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {e}")
        exit(1)


def align_reads(ref_genome, fastq_file1, fastq_file2, output_sam):
    """ Align reads to the reference genome """
    print("Aligning reads to the reference genome...")

    bwa_command = ['bwa', 'mem', ref_genome, fastq_file1, fastq_file2]
    run_command(bwa_command, output_sam)


def convert_sam_to_bam(sam_file, bam_file):
    """ Convert SAM file to BAM file """
    print("Converting SAM to BAM...")
    run_command(['samtools', 'view', '-b', sam_file, '-o', bam_file])

def sort_bam(bam_file, sorted_bam_file):
    """ Sort BAM file """
    print("Sorting BAM file...")
    run_command(['samtools', 'sort', bam_file, '-o', sorted_bam_file])

def index_bam(sorted_bam_file):
    """ Index the sorted BAM file """
    print("Indexing the sorted BAM file...")
    run_command(['samtools', 'index', sorted_bam_file])

def main():
    parser = argparse.ArgumentParser(description='Process sequencing reads.')
    parser.add_argument('--ref_genome', required=True, help='Reference genome file (FASTA)')
    parser.add_argument('--fastq1', required=True, help='First FASTQ file with sequencing reads')
    parser.add_argument('--fastq2', required=True, help='Second FASTQ file with sequencing reads')
    parser.add_argument('--output', required=True, help='Prefix for output files')

    args = parser.parse_args()

    sam_file = f"{args.output}.sam"
    bam_file = f"{args.output}.bam"
    sorted_bam_file = f"{args.output}_sorted.bam"

    align_reads(args.ref_genome, args.fastq1, args.fastq2, sam_file)
    convert_sam_to_bam(sam_file, bam_file)
    sort_bam(bam_file, sorted_bam_file)
    index_bam(sorted_bam_file)

if __name__ == "__main__":
    main()
