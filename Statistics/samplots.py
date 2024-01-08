# Must be run using Python 2.7.x


import subprocess
import os

def run_samplot(vcf_path, bam_path, output_dir, tag):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue  # Skip header lines

            columns = line.strip().split('\t')
            chrom, pos, sv_id, ref, alt, _, _, info = columns[:8]
            info_dict = dict((k, v) for k, v in (field.split('=') for field in info.split(';') if '=' in field))

            sv_type = info_dict.get('SVTYPE', 'UNK')
            end = info_dict.get('END', pos)

            # Constructing samplot command
            output_file = "{}/{}_{}_{}_{}_{}.png".format(output_dir, tag, sv_id, chrom, pos, end)
            cmd = [
                "samplot", "plot",
                "-n", sv_id,
                "-b", bam_path,
                "-o", output_file,
                "-c", chrom,
                "-s", pos,
                "-e", end,
                "-t", sv_type
            ]

            # Run samplot command
            print "Running: {}".format(' '.join(cmd))
            subprocess.check_call(cmd)

# Define paths
candidate_vcf_path = "/home/luke/Desktop/TargetMI/SV_Analysis/smoove_output/N12878_candidates-smoove.genotyped.vcf"
# somatic_vcf_path = "/home/e/Documents/manta/SV_finder_tools/results/Data/somaticSV.vcf"
bam_path = "/home/luke/Desktop/TargetMI/SV_Analysis/bam_file/N12878_bwamem_bam.bam"
output_dir = "/home/luke/Desktop/TargetMI/SV_Analysis/samplot"

# Run for both VCF files
run_samplot(candidate_vcf_path, bam_path, output_dir, "candidate")
# run_samplot(somatic_vcf_path, bam_path, output_dir, "somatic")