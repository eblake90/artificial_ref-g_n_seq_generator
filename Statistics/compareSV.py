import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
import os
import argparse

# Setting up Argument Parser
parser = argparse.ArgumentParser(description='Process VCF and BAM files for SV analysis.')
parser.add_argument('--bam', help='Path to the BAM file', default=None)
parser.add_argument('--vcf', required=True, help='Path to the VCF file')
parser.add_argument('--output', required=True, help='Path to the output folder')
args = parser.parse_args()

# File paths
bam_file = args.bam
vcf_file = args.vcf
output_folder = args.output

# Extract the file name from the VCF file path
vcf_file_name = os.path.basename(vcf_file).split('.')[0]
bam_file_name = os.path.basename(bam_file).split('.')[0] if bam_file else "NoBAM"

# Create Data and Graphs folders in the specified output directory
save_graph_path = os.path.join(output_folder, 'Graphs')
save_data_path = os.path.join(output_folder, 'Data')
os.makedirs(save_graph_path, exist_ok=True)
os.makedirs(save_data_path, exist_ok=True)

# Function to get read depth from a BAM file
def get_read_depth(bam_path):
    read_depth = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileupcolumn in bam.pileup():
            read_depth.append(pileupcolumn.n)
    return read_depth

# Function to categorize SVs
def categorize_sv(record):
    if 'SVTYPE' in record.info:
        return record.info['SVTYPE']
    return 'Unknown'

# Count SV types for VCF
def count_sv_types(vcf_path):
    sv_counts = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            sv_type = categorize_sv(record)
            sv_counts[sv_type] = sv_counts.get(sv_type, 0) + 1
    return sv_counts

# Count SVs in VCF file
sv_counts = count_sv_types(vcf_file)

# Read depth related variables
read_depth_data = []
average_read_depth = lowest_read_depth = highest_read_depth = None

# Get read depth data only if BAM file is provided
if bam_file:
    read_depth_data = get_read_depth(bam_file)
    average_read_depth = sum(read_depth_data) / len(read_depth_data)
    lowest_read_depth = min(read_depth_data)
    highest_read_depth = max(read_depth_data)

# Convert data to pandas DataFrame for easy handling
sv_counts_df = pd.DataFrame(list(sv_counts.items()), columns=['SV Type', 'Count'])
read_depth_df = pd.DataFrame(read_depth_data, columns=['Read Depth']) if read_depth_data else pd.DataFrame()

# Save DataFrames to CSV
sv_counts_csv_path = os.path.join(save_data_path, f'{vcf_file_name}_sv_counts.csv')
sv_counts_df.to_csv(sv_counts_csv_path, index=False)
if read_depth_data:
    read_depth_csv_path = os.path.join(save_data_path, f'{bam_file_name}_read_depth.csv')
    read_depth_df.to_csv(read_depth_csv_path, index=False)

# Plotting SV types for VCF
plt.figure(figsize=(10, 5))
sns.barplot(x=list(sv_counts.keys()), y=list(sv_counts.values()))
plt.xlabel('SV Type')
plt.ylabel('Count')
plt.title(f'Counts of Structural Variant Types from {vcf_file_name}')
plt.xticks(rotation=45)
plot_file_path = os.path.join(save_graph_path, f'{vcf_file_name}_sv_counts.png')
plt.savefig(plot_file_path)

# Plotting read depth if BAM file is provided
if read_depth_data:
    plt.figure(figsize=(12, 6))
    plt.plot(read_depth_data, color='blue', linewidth=0.5)
    plt.xlabel('Genomic Position')
    plt.ylabel('Read Depth')
    plt.title(f'Read Depth Distribution of {bam_file_name}')
    read_depth_plot_path = os.path.join(save_graph_path, f'{bam_file_name}_read_depth.png')
    plt.savefig(read_depth_plot_path)

# Convert SV counts DataFrame to a list of lists for Table
sv_counts_table_data = [['SV Type', 'Count']] + sv_counts_df.values.tolist()
# Constants for image sizing and positioning
image_width = 400
image_height = 200
margin = 50
# Creating a PDF file
pdf_file_path = os.path.join(output_folder, f"{vcf_file_name}_summary_report.pdf")
pdf = canvas.Canvas(pdf_file_path, pagesize=letter)
width, height = letter

# Title and Subtitle
pdf.setFont("Helvetica-Bold", 16)
pdf.drawString(100, height - 40, vcf_file_name)
pdf.setFont("Helvetica", 12)
pdf.drawString(100, height - 60, bam_file_name)

# Adding text and images to PDF
if os.path.exists(plot_file_path):
    pdf.drawImage(plot_file_path, 100, height - margin - image_height, width=image_width, height=image_height)
else:
    print(f"Warning: Image file not found at {plot_file_path}")

pdf.showPage()

if os.path.exists(read_depth_plot_path):
    pdf.drawImage(read_depth_plot_path, 100, height - margin - image_height, width=image_width, height=image_height)
else:
    print(f"Warning: Image file not found at {read_depth_plot_path}")

pdf.showPage()  # Start a new page for tables

# Adding Read Depth Statistics table
pdf.setFont("Helvetica-Bold", 14)
pdf.setFillColor(colors.black)  # Explicitly set text color to black
pdf.drawString(100, height - 100, "Read Depth Statistics")
read_depth_stats_table = Table([
    ['Metric', 'Value'],
    ['Average', round(average_read_depth, 2)],
    ['Lowest', lowest_read_depth],
    ['Highest', highest_read_depth]
], colWidths=[200, 100])
read_depth_stats_table.setStyle(TableStyle([
    ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
    ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
    ('BOX', (0,0), (-1,-1), 1, colors.black),
    ('INNERGRID', (0,0), (-1,-1), 1, colors.black),
]))
read_depth_stats_table.wrapOn(pdf, width, height)
read_depth_stats_table.drawOn(pdf, 100, height - 150)

# Adding SV Counts table
pdf.setFont("Helvetica-Bold", 14)
pdf.drawString(100, height - 300, "SV Counts")
sv_counts_table = Table(sv_counts_table_data, colWidths=[200, 100])
sv_counts_table.setStyle(TableStyle([
    ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
    ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
    ('BOX', (0,0), (-1,-1), 1, colors.black),
    ('INNERGRID', (0,0), (-1,-1), 1, colors.black),
]))
sv_counts_table.wrapOn(pdf, width, height)
sv_counts_table.drawOn(pdf, 100, height - 350)

# Save PDF
pdf.save()

print(f"PDF report saved at: {pdf_file_path}")