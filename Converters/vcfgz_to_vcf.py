import os
import gzip
import shutil
import argparse

def unzip_and_move_vcf_files(source_dir, destination_dir):
    if not os.path.exists(source_dir):
        print("Source directory '{}' does not exist.".format(source_dir))
        exit(1)

    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    for root, _, files in os.walk(source_dir):
        for filename in files:
            if filename.endswith(".vcf.gz"):
                vcf_gz_file = os.path.join(root, filename)
                vcf_filename, _ = os.path.splitext(filename)
                vcf_filename, _ = os.path.splitext(vcf_filename)
                with gzip.open(vcf_gz_file, 'rb') as f_in, open(
                        os.path.join(destination_dir, "{}.vcf".format(vcf_filename)), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                print("Extracted and moved '{}.vcf'".format(vcf_filename))

    print("Finished.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Unzip and move VCF files.')
    parser.add_argument('--input', required=True, help='Path to the folder containing .vcf.gz files')
    parser.add_argument('--output', required=True, help='Path to the output folder for .vcf files')
    args = parser.parse_args()

    unzip_and_move_vcf_files(args.input, args.output)
