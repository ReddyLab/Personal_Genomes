#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:37:20 2024

@author: kd259
"""

# %% Import packages

import pysam
import argparse
import gzip

# %% Define inputs
parser = argparse.ArgumentParser(
                    prog='count_alleles',
                    description='This program counts the reference and alternate alleles in each fragment',
                    epilog='')
parser.add_argument("-s", "--sam_name", help="The path to the sam file to search for variants")
parser.add_argument("-o", "--output", help="The path to the processed file aligned to GRCh38")
parser.add_argument("-v", "--vcf_name", help="The path to the VCF file containing the variants in this sam file")
parser.add_argument("-r", "--inread", action='store_true', help="Flag to indicate that the variants should be checkin only in the read instead of the whole fragment")

args = parser.parse_args()

# Define functions
# this function determines the start end end coordinate, as well as the alt alleles
# It takes mate_names, and the sam file line
# This also needs to work with pysam
# The output is going to be the start coordinate, end coordiante, and List
# I need to make sure that I have soemthing in there for insertions that start before the fragment
#   Specifically, I need to make sure that they don't get deleted from the vcf list file too early
# Also, I need to make sure that the whole last chromosome is printed, because that's an issue sometimes
# Actually, just make sure that the end of every chromosome is printed
# Also, I need to test this on data that has been manually typed
def find_end(start, vars, length):
    end_val = start + length
    for v in vars:
        end_val += len(v[3]) - len(v[2])
    return end_val

if args.inread:
    # This is the function for if you only use the variants in the reads themselves
    def find_variants(sam_line, seen_fragment_names):
        # Check if the other read in this fragment has been used
        if sam_line.qname in seen_fragment_names:
            # If the other read in this fragment has been used, do not use this one
            mate_ind = seen_fragment_names.index(sam_line.qname)
            delete = seen_fragment_names.pop(mate_ind)
            return "SKIP"
        else:
            # Otherwise, figure out which method to use to include this read
            # Determine overlap
            start_pos = sam_line.pos
            m_start_pos = sam_line.mpos
            if sam_line.has_tag("rv"): vrs = sam_line.get_tag("rv").split(";") # If there are fragment variants, retrieve them
            else: vrs = [] # Otherwise, create a blank list
            end_pos = find_end(start_pos, vrs, abs(sam_line.rlen))
            # Case 1: no overlap between reads (use just the read variants and do the same for the other read) (first time encountering the fragment)
            if end_pos < m_start_pos and start_pos < m_start_pos:
                return start_pos, end_pos, vrs # Return the values
            # Case 2: some or all overlap between reads (use the fragment variants and do not use the other read)
            elif start_pos <= m_start_pos and end_pos >= m_start_pos:
                if sam_line.is_paired: seen_fragment_names.append(sam_line.qname) # Save the read name if it is paired
                if sam_line.has_tag("fv"): vrs = sam_line.get_tag("fv").split(";") # If there are fragment variants, retrieve them
                else: vrs = [] # Otherwise, create a blank list
                start_val = sam_line.pos # Recover the start coordinate
                end_val = find_end(start_val, vrs, abs(sam_line.tlen)) # Recover the end coordinate
                return start_val, end_val, vrs # Return the values
            # Case 1: no overlap between reads (second time encountering the fragment)
            elif m_start_pos < start_pos:
                return start_pos, end_pos, vrs # Return the values

else:
    # This is the function for if you use the variants in the whole fragment
    def find_variants(sam_line, seen_fragment_names):
        # Check if the other read for this fragment has been used
        if sam_line.qname in seen_fragment_names:
            # If the other read for this fragment has been used, do not use this one
            mate_ind = seen_fragment_names.index(sam_line.qname)
            delete = seen_fragment_names.pop(mate_ind)
            return "SKIP"
        else:
            # Otherwise, find the start and end coordinates, and the variants in the whole fragment
            if sam_line.is_paired: seen_fragment_names.append(sam_line.qname) # Save the read name if it is paired
            if sam_line.has_tag("fv"): vrs = sam_line.get_tag("fv").split(";") # If there are fragment variants, retrieve them
            else: vrs = [] # Otherwise, create a blank list
            start_val = sam_line.pos # Recover the start coordinate
            end_val = find_end(start_val, vrs, abs(sam_line.tlen)) # Recover the end coordinate
            return start_val, end_val, vrs # Return the values

# Load VCF
vcf_file = gzip.open(args.vcf_name, 'rb')
vcf = []
for i in range(23):
    vcf.append([])

for line in vcf_file:
    line_split = str(line)[2:-3].replace("\\t", "\t").split("\t")
    if "X" in line_split[0]: chrom = 23
    else: chrom = int(line_split[0][3:])
    important_part = [line_split[0], line_split[1], line_split[3], line_split[4], 0, 0]
    vcf[chrom-1].append(important_part)

# Check sam file
sam_file = pysam.AlignmentFile(args.sam_name, "rb")
chrom = 0 # Index of the current chromosome
variant = 0 # Index of the variant in the chromosome
vcf_write = [] # List of variants to write to the new vcf file
new_vcf = open(args.output, 'w')
mate_names = []
for line in sam_file:

    chromosome = line.reference_name
    # If the chromosome changes
    if chromosome[3:] not in ["M", "X", "Y"] and int(chromosome[3:])-1 != chrom:
        # Write everything to the write vcf with the counts
        while variant < len(vcf[chrom]):
            vcf_write.append(vcf[chrom][variant])
            variant += 1
        while len(vcf_write) > 0:
            vcf_line = vcf_write.pop(0)
            vcf_line[-2] = str(int(vcf_line[-2]))
            vcf_line[-1] = str(int(vcf_line[-1]))
            vcf_line.append(str(int(vcf_line[-2]) + int(vcf_line[-1])))
            new_vcf.write("\t".join(vcf_line) + "\n")
        # reset the chromosome and the position index values
        chrom = int(chromosome[3:])-1
        variant = 0
        mate_names = []

    # If the chromosome changes to the X chromosome
    elif chromosome[3:] in "X" and chrom != 22:
        # Write everything to the write vcf with the counts
        while variant < len(vcf[chrom]):
            vcf_write.append(vcf[chrom][variant])
            variant += 1
        while len(vcf_write) > 0:
            vcf_line = vcf_write.pop(0)
            vcf_line[-2] = str(int(vcf_line[-2]))
            vcf_line[-1] = str(int(vcf_line[-1]))
            vcf_line.append(str(int(vcf_line[-2]) + int(vcf_line[-1])))
            new_vcf.write("\t".join(vcf_line) + "\n")
        # reset the chromosome and the position index values
        chrom = 22
        variant = 0
        mate_names = []

    # Check for if you've moved past the first thing
    check_pos = 0
    while len(vcf_write) > check_pos and int(vcf_write[check_pos][1]) + len(vcf_write[check_pos][2]) - 1 < int(line.pos):
        vcf_line = vcf_write.pop(0)
        vcf_line[-2] = str(int(vcf_line[-2]))
        vcf_line[-1] = str(int(vcf_line[-1]))
        vcf_line.append(str(int(vcf_line[-2]) + int(vcf_line[-1])))
        new_vcf.write("\t".join(vcf_line) + "\n")

    vcf_info = find_variants(line, mate_names)
    if vcf_info != "SKIP":
        # Separate the information for the vcf
        start_val, end_val, vrs = vcf_info

        while variant < len(vcf[chrom]) and int(vcf[chrom][variant][1]) < end_val:
            vcf_write.append(vcf[chrom][variant])
            variant += 1

        pos_s = []
        for v in vrs:
            pos_s.append(v[1])

        # Add counts to the new vcf
        adjustment = 0
        for v in vcf_write:
            if int(v[1]) + adjustment < end_val:
                v_string = v[0] + ":" + v[1] + "," + v[2] + "," + v[3]
                pos = v[1]
                if v_string in vrs:
                    v[-1] += 1
                    adjustment += len(v[3]) - len(v[2])
                elif int(v[1]) + len(v[2]) > start_val and pos not in pos_s:
                    v[-2] += 1
                    #adjustment += min(len(v[3]) - len(v[2]), 0)

while variant < len(vcf[chrom]):
    vcf_write.append(vcf[chrom][variant])
    variant += 1


for line in vcf_write:
    line[-2] = str(int(line[-2]))
    line[-1] = str(int(line[-1]))
    line.append(str(int(line[-2]) + int(line[-1])))
    new_vcf.write("\t".join(line) + "\n")
