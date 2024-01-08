# run in py3
# to run
# python fasta_to_fa.py input.fasta output.fa


import argparse


def convert_fasta_to_fa(input_file, output_file):
    # Open the input FASTA file for reading.
    with open(input_file, 'r') as fasta_file:
        # Read the content of the file.
        data = fasta_file.read()

    # Open the output FA file for writing.
    with open(output_file, 'w') as fa_file:
        # Write the content to the FA file.
        fa_file.write(data)


def main():
    # Set up command line argument parsing.
    parser = argparse.ArgumentParser(description='Convert FASTA file to FA file.')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output FA file')

    # Parse the command line arguments.
    args = parser.parse_args()

    # Call the conversion function.
    convert_fasta_to_fa(args.input, args.output)


if __name__ == "__main__":
    main()
