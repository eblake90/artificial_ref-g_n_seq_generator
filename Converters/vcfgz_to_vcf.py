import os
import gzip
import shutil

# Source and destination directories
source_dir = "/home/eblake/Documents/manta/SV_finder_tools/MantaWorkflow/results/variants"
destination_dir = "/home/eblake/Documents/manta/SV_finder_tools/results/Data"

# Check if the source directory exists
if not os.path.exists(source_dir):
    print("Source directory '{}' does not exist.".format(source_dir))
    exit(1)

# Check if the destination directory exists; create it if it doesn't
if not os.path.exists(destination_dir):
    os.makedirs(destination_dir)

# Loop through .vcf.gz files in the source directory
for root, _, files in os.walk(source_dir):
    for filename in files:
        if filename.endswith(".vcf.gz"):
            vcf_gz_file = os.path.join(root, filename)
            # Extract the filename without the path and double extension
            vcf_filename, _ = os.path.splitext(filename)
            vcf_filename, _ = os.path.splitext(vcf_filename)

            # Unzip the .vcf.gz file to .vcf format and move it to the destination directory
            with gzip.open(vcf_gz_file, 'rb') as f_in, open(
                    os.path.join(destination_dir, "{}.vcf".format(vcf_filename)), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            print("Extracted and moved '{}.vcf'".format(vcf_filename))

print("Finished.")