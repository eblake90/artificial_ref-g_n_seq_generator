# to run
# python vcf_to_csv_converter.py path/to/input.vcf path/to/output.csv

import csv
import argparse

def convert_vcf_to_csv(vcf_file_path, csv_file_path):
    with open(vcf_file_path, 'r') as vcf_file, open(csv_file_path, 'w', newline='') as csv_file:
        vcf_reader = csv.reader(vcf_file, delimiter='\t')
        csv_writer = csv.writer(csv_file)

        for row in vcf_reader:
            if row[0].startswith('##'):
                continue
            elif row[0].startswith('#'):
                # Remove '#' from the first column header
                row[0] = row[0][1:]
                csv_writer.writerow(row)
            else:
                csv_writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description='Convert a VCF file to a CSV file.')
    parser.add_argument('vcf_file', type=str, help='Path to the input VCF file')
    parser.add_argument('csv_file', type=str, help='Path to the output CSV file')

    args = parser.parse_args()
    convert_vcf_to_csv(args.vcf_file, args.csv_file)

if __name__ == '__main__':
    main()
