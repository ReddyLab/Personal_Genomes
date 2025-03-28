#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import packages
import sys
import gzip
from Pipe import Pipe
import math
import numpy as np

# %% Define inputs
parser = argparse.ArgumentParser(
                    prog='write_genome_and_all_vcf',
                    description='This program creates a personal genome and writes the corresponding VCFs',
                    epilog='Text at the bottom of help')
parser.add_argument("-s", "--samples", help="A comma separated list of sample IDs")
parser.add_argument("-c", "--chromosome", help="The chromosome to create")
parser.add_argument("-f", "--frag_len", default=1000, type=int, help="The maximum length of a fragment")
parser.add_argument("-v", "--vcf_name", help="The path to the gzipped VCF file")
parser.add_argument("--twoBitToFa", help="The path to the twoBitToFa program")
parser.add_argument("--twobit", help="The path to the twobit reference genome file")

parser.add_argument("-o", "--output", help="The path to the processed file aligned to GRCh38")

args = parser.parse_args()


samples_string = args.samples # "NA18517,NA18519,NA18858,NA18868,NA18873"
samples_list = samples_string.split(",")
chrom_in = args.chromosome
fasta_filename = "chr" + chrom_in + ".fasta"
temp_filename = "chr" + chrom_in + ".temp"
fragment_length = args.frag_len
vcf_name = args.vcf_name

# %% Define the chromosome lengths

chr_lengths = [248956422, # chr1
               242193529, # chr2
               198295559, # chr3
               190214555, # chr4
               181538259, # chr5
               170805979, # chr6
               159345973, # chr7
               145138636, # chr8
               138394717, # chr9
               133797422, # chr10
               135086622, # chr11
               133275309, # chr12
               114364328, # chr13
               107043718, # chr14
               101991189, # chr15
               90338345,  # chr16
               83257441,  # chr17
               80373285,  # chr18
               58617616,  # chr19
               64444167,  # chr20
               46709983,  # chr21
               50818468,  # chr22
               156040895  # chrX
               ]


# %% Define functions
def process_indel(indel_list, vcf_file, het_vcf, fasta_file, add_distance):

    # Define the indel itself and all variants in LD
    indel_row = indel_list[0]
    indel_close = indel_list[1]
    temp_filename = "temp." + indel_row[0] + ".txt"

    heterozygous = False
    # Check that the variant is heterozygous
    if not (indel_row[5:] == ["0|0"] * (len(indel_row)-5) or indel_row[5:] == ["1|1"] * (len(indel_row)-5)):
        heterozygous = True
        het_vcf.write("\t".join(indel_row) + "\n")

        # Check that the variants in LD are homozygous
        all_homozygous = True
        alt_indel = False
        het_indel = False
        alt_variants = []
        for row in indel_close:
            if not (row[5:] == ["0|0"] * (len(row)-5) or row[5:] == ["1|1"] * (len(row)-5)):
                all_homozygous = False
                if len(row[3]) != len(row[4]):
                    het_indel = True
                else:
                    pass
                    # This is where I should phase the variants to the alt or ref allele of the thing and then see which sequence each thing is associated with
            if (row[5:] == ["1|1"] * (len(row)-5)) and (len(row[3]) != len(row[4])):
                alt_indel = True
            if row[5:] == ["1|1"] * (len(row)-5):
                alt_variants.append(row)

    # Find the hg38 reference sequence
    if heterozygous and all_homozygous and not alt_indel and not het_indel:
        chrom = indel_row[0]
        start_pos = str(int(indel_row[1])-add_distance)
        end_pos = str(int(indel_row[1])+add_distance + 1)
        ref = Pipe.run(args.twoBitToFa + " -seq="+chrom+" -start="+start_pos+" -end="+end_pos+" " + args.twobit + " "+temp_filename)

        # Load the new file
        temp_file = open(temp_filename)
        sequence_name = temp_file.readline()[:-1]
        sequence = ""
        for line in temp_file:
            sequence += line[:-1]
        ref_sequence = sequence
        # Define the refernce sequence for these variants
        for var in alt_variants:
            if len(var[3]) == 1 and len(var[4]) == 1:
                variant_pos = abs(int(var[1]) - int(start_pos)) - 1
                ref_sequence = ref_sequence[:variant_pos] + var[4] + ref_sequence[variant_pos+1:]

        # Save the reference sequence for these variants
        fasta_file.write(sequence_name + "_ref\n")
        fasta_file.write(ref_sequence + "\n")

        # Define the alternate sequence for these variants
        variant_pos = 249
        alt_sequence = ref_sequence[:variant_pos] + indel_row[4] + ref_sequence[variant_pos+len(indel_row[3]):]

        # Save the alternate sequence for these variants
        fasta_file.write(sequence_name + "_alt\n")
        fasta_file.write(alt_sequence + "\n")

        # Save to the overall VCF file
        vcf_file.write("\t".join(indel_row) + "\n")



def process_new_sequence(fasta_file, variants, chrom, start_pos, end_pos, temp_filename):

    # Find the hg38 reference sequence
    ref = Pipe.run(args.twoBitToFa + " -seq=chr"+chrom+" -start="+start_pos+" -end="+end_pos+" " + args.twobit + " "+temp_filename)

    # Load the new file to get the sequence
    temp_file = open(temp_filename)
    source_sequence_name = temp_file.readline()[:-1]
    ref_sequence = ""
    for line in temp_file:
        ref_sequence += line[:-1]


    # Define the different options for the haplotypes of this set of variants
    if len(variants) > 0:
        all_haplotypes = []
        for ind in range(len(variants[0])-9):
            new_hap = ""
            for v in range(len(variants)):
                new_hap += variants[v][ind+9][0]
            if new_hap not in all_haplotypes:
                all_haplotypes.append(new_hap)
            new_hap = ""
            for v in range(len(variants)):
                if len(variants[v][ind+9]) > 1:
                    new_hap += variants[v][ind+9][2]
                else:
                    new_hap += variants[v][ind+9][0]
            if new_hap not in all_haplotypes:
                all_haplotypes.append(new_hap)


        # Edit the sequence once for each haplotype in this set
        for hap_str in all_haplotypes:

            # Initialize a new header and a new sequence
            new_head = source_sequence_name
            new_seq = ref_sequence

            # Edit the header for all variants to be included
            for i in range(len(hap_str)):
                if hap_str[i] == "1":
                    new_head += ";" + chrom + ":" + variants[i][1]
                    new_head += "," + variants[i][3] + "," + variants[i][4]

            # Change all the indels in the sequence in reverse order
            for i in range(len(hap_str)):
                i_r = (i + 1) * -1
                if hap_str[i_r] == "1":
                    variant_pos = abs(int(variants[i_r][1]) - int(start_pos)) - 1
                    new_seq = new_seq[:variant_pos] + variants[i_r][4] + new_seq[variant_pos+len(variants[i_r][3]):]

            # Save the sequence
            fasta_file.write(new_head + "\n")
            fasta_file.write(new_seq + "\n")

    else:
        # Initialize a new header and a new sequence
        new_head = source_sequence_name
        new_seq = ref_sequence
        # Save the sequence
        fasta_file.write(new_head + "\n")
        fasta_file.write(new_seq + "\n")



# %% Create header
vcf_file = gzip.open(vcf_name)

# Find the header line
target_chr = "##" # This character indicates a header line
next_line = "" # Initialize the next line
while target_chr == "##": # While the previou line was a header line
    prev_line = next_line # Save the previous line
    next_line = vcf_file.readline() # Read the next line in the file
    target_chr = str(next_line)[2:4] # Save the target character in this new line
# Delete the target character because it won't be used again
del target_chr


# %% Define the important columns in the file

# Find the locations of the important information in each line
# The first 5 columns have general info about the SNP that are important
# Then you need the column numbers for each sample in the populations that were specified in the input

# Save the header of the file as a list
header_list = str(next_line).split("\\t") # Split the header line into a list
total_samples = len(header_list[9:])

important_locations = [0,1,2,3,4,5,6,7,8] # Save the chromosome and position, and the alt and ref sequences
condensed_header_locations = [0,1,2,3,4,5,6,7,8]

for sample in samples_list: # For each sample in the pool
    try:
        # Put the column number of that sample into the list of columns if it is in the file
        important_locations.append(header_list.index(sample))
        # Save a shorter list with just the samples included in the important locations list
        condensed_header_locations.append(sample)
    except:
        pass



# %% Search the VCF file

vcf_list = []
all_indels = []
unfinished_indels = []
temp_list = []
homozygous_alt = []

all_indels_vcf_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
all_indels_vcf_header += samples_list

vcf_out = open("VCF/chr" + chrom_in + ".vcf", 'w')
personal_vcf_files, personal_vcf_indices = [], []
for f in range(len(samples_list)):
    personal_vcf_files.append(open(samples_list[f] + "/chr" + chrom_in + ".vcf", 'w'))
    personal_vcf_indices.append(important_locations[all_indels_vcf_header.index(samples_list[f])])
fasta_file = open("Genome/chr" + chrom_in + ".fasta", 'w')
next_start = "1"
First = True

for line in vcf_file:
  if len(str(line)) > 3:
    line_list = str(line)[2:-3].split("\\t")
    small_row = []
    for ind in important_locations:
        small_row.append(line_list[ind])

    # Write to the general vcf
    if small_row[9:] == ["1|1"] * (len(samples_list)):
        vcf_out.write("\t".join(small_row) + "\n")
    elif not small_row[9:] == ["0|0"] * (len(samples_list)):
        vcf_out.write("\t".join(small_row) + "\n")

    # Write to the fasta file
    #Add alt alleles to the list of alt variants
    if not small_row[9:] == ["0|0"] * (len(small_row)-9):
        temp_list.append(small_row)

    # Check for max variants fragment ending
    if len(temp_list) > max_vars and (int(temp_list[-1][1]) - 1 - int(next_start)) > (2*fragment_length):
        # Find the end coordinate
        end_coord = str( int(temp_list[-1][1]) -1 )
        # Paste the sequence in
        process_new_sequence(fasta_file, temp_list[:-1], chrom_in, next_start, end_coord, temp_filename)
        # Re-initialize the start coordinate and the list of variants
        next_start = str(int(end_coord) - fragment_length - 1)
        temp_list.reverse()
        temp_list_new = []
        i = 0
        while int(temp_list[i][1]) > int(next_start):
            temp_list_new.append(temp_list[i])
            i += 1
        temp_list_new.reverse()
        temp_list = temp_list_new
        del temp_list_new

    # Check for constant region fragment ending
    elif len(temp_list) > 1 and (int(temp_list[-1][1]) - int(temp_list[-2][1])) > fragment_length + 2:
        # Process the first variable region
        # find the end coordinate
        end_coord = str( int(temp_list[-2][1]) + fragment_length + 1)
        # Paste the sequence in
        process_new_sequence(fasta_file, temp_list[:-1], chrom_in, next_start, end_coord, temp_filename)
        # Re-initialize the start coordinate
        next_start = str( int(temp_list[-2][1]) + 1 )
        # Process the constant region
        # find the end coordinate
        end_coord = str( int(temp_list[-1][1]) - 1)
        # Paste the sequence in
        process_new_sequence(fasta_file, [], chrom_in, next_start, end_coord, temp_filename)
        # Re-initialize the start coordinate and the list of variants
        next_start = str(int(end_coord) - fragment_length - 1)
        temp_list = [temp_list[-1]]

    # Check for telomere fragment ending
    elif First and len(temp_list) > 0:
        First = False
        # Find the end coordinate
        end_coord = str( int(temp_list[0][1]) -1 )
        # Paste the sequence in
        process_new_sequence(fasta_file, [], chrom_in, next_start, end_coord, temp_filename)
        # Re-initialize the start coordinate and the list of variants
        next_start = str(int(end_coord) - fragment_length - 1)

    # Check for a personal variant
    personal_variant_check = np.unique(line_list[9:], return_counts=True)
    if len(personal_variant_check[0]) == 2:
        if "0|0" == personal_variant_check[0][0] and personal_variant_check[1][0] == total_samples-1 and not small_row[9:] == ["0|0"] * (len(small_row)-9):
            vcf_file_index = small_row[9:].index(personal_variant_check[0][1])
            personal_vcf_files[vcf_file_index].write("\t".join(small_row) + "\n")
        if "0|0" == personal_variant_check[0][1] and personal_variant_check[1][1] == total_samples-1 and not small_row[9:] == ["0|0"] * (len(small_row)-9):
            vcf_file_index = small_row[9:].index(personal_variant_check[0][0])
            personal_vcf_files[vcf_file_index].write("\t".join(small_row) + "\n")
        if "1|1" == personal_variant_check[0][0] and personal_variant_check[1][0] == total_samples-1 and not small_row[9:] == ["1|1"] * (len(small_row)-9):
            vcf_file_index = small_row[9:].index(personal_variant_check[0][1])
            personal_vcf_files[vcf_file_index].write("\t".join(small_row) + "\n")
        if "1|1" == personal_variant_check[0][1] and personal_variant_check[1][1] == total_samples-1 and not small_row[9:] == ["1|1"] * (len(small_row)-9):
            vcf_file_index = small_row[9:].index(personal_variant_check[0][1])
            personal_vcf_files[vcf_file_index].write("\t".join(small_row) + "\n")

# Save the vcf file
vcf_out.close()

# Save the personal vcf files
for f in personal_vcf_files:
    f.close()

# Write the last set to the fasta file
end_coord = str( int(temp_list[-1][1]) + fragment_length + 1)
try:
    chrom_end = str(chr_lengths[int(chrom_in)-1])
except:
    chrom_end = str(chr_lengths[-1])

if int(end_coord) < int(chrom_end) - fragment_length:
    end_coord = str(int(temp_list[-1][1]) + 1)
    process_new_sequence(fasta_file, temp_list, chrom_in, next_start, end_coord, temp_filename)
    next_start = str( int(end_coord) - fragment_length - 1)
    temp_list = []

process_new_sequence(fasta_file, temp_list, chrom_in, next_start, chrom_end, temp_filename)
