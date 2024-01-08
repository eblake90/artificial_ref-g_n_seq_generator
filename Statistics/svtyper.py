# Must be run using Python 2.7.x

import svtyper.classic as svt

input_vcf = "/home/luke/Desktop/TargetMI/SV_Analysis/smoove_output/N12878_candidates-smoove.genotyped.vcf"
input_bam = "/home/luke/Desktop/TargetMI/SV_Analysis/bam_file/N12878_bwamem_bam.bam"
library_info = "/home/luke/Desktop/TargetMI/SV_Analysis/bam_file/N12878_bwamem_bam.bam.json"
output_vcf = "/home/luke/Desktop/TargetMI/SV_Analysis/smoove_output/N12878_candidates-genotyped.vcf"

with open(input_vcf, "r") as inf, open(output_vcf, "w") as outf:
    svt.sv_genotype(bam_string=input_bam,
                    vcf_in=inf,
                    vcf_out=outf,
                    min_aligned=20,
                    split_weight=1,
                    disc_weight=1,
                    num_samp=1000000,
                    lib_info_path=library_info,
                    debug=False,
                    alignment_outpath=None,
                    ref_fasta=None,
                    sum_quals=False,
                    max_reads=None,
                    max_ci_dist=None)