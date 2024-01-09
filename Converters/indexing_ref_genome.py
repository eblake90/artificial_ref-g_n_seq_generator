# run on py3
# to run
# python /home/eblake/Documents/github/1.1/artificial_ref-g_n_seq_generator/Converters/indexing_ref_genome.py --ref_genome ref_genome.fa

import subprocess
import argparse

def run_command(command):
    """ Run a shell command """
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {e}")
        exit(1)

def index_genome(ref_genome):
    """ Index the reference genome using BWA """
    print("Indexing the reference genome...")
    run_command(['bwa', 'index', ref_genome])

def create_fai_file(ref_genome):
    """ Create a .fa.fai index file for the reference genome """
    print(f"Creating .fa.fai index file for {ref_genome}...")
    run_command(['samtools', 'faidx', ref_genome])

def main():
    parser = argparse.ArgumentParser(description='Index a reference genome and create .fa.fai file.')
    parser.add_argument('--ref_genome', required=True, help='Reference genome file (FASTA)')

    args = parser.parse_args()

    index_genome(args.ref_genome)
    create_fai_file(args.ref_genome)

if __name__ == "__main__":
    main()
