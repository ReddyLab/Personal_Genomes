#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 12:56:20 2023

@author: kd259
"""

# %% Import packages

import csv
import argparse
import file_headers

# %% Define inputs
parser = argparse.ArgumentParser(
                    prog='convert_to_hg38',
                    description='This program converts the personal genome alignment to an hg38 alignment',
                    epilog='Text at the bottom of help')
parser.add_argument("-i", "--input", help="The path to the unprocessed file aligned to the personal genome")
parser.add_argument("-o", "--output", help="The path to the processed file aligned to GRCh38")

args = parser.parse_args()


# %% Define functions

def find_name_and_vars(in_string, start_ind, frag_len):

    in_list = in_string.split(";")

    seq_name = in_list[0]
    variants = tuple(in_list[1:])
    frag_start = int(seq_name[seq_name.index(":")+1:seq_name.index("-")])

    var_locs = []
    adjustment = 0
    variants_short, var_locs_short = [], []
    for var in variants:
        var_list = var.split(",")
        position = int(var[var.index(":")+1:var.index(",")])

        if len(var_locs) > 0:
            new_loc = max(position - frag_start + 1 + adjustment, var_locs[-1])
        else:
            new_loc = position - frag_start + 1 + adjustment

        adjustment += len(var_list[2]) - len(var_list[1])

        var_locs.append(new_loc)

        if (new_loc + len(var_list[2])) > start_ind and new_loc < start_ind + frag_len:
            variants_short.append(var)
            var_locs_short.append(new_loc)


    var_locs = tuple(var_locs)

    variants_short = tuple(variants_short)
    var_locs_short = tuple(var_locs_short)

    return seq_name, variants_short, var_locs_short, variants, var_locs

def process_variants(vcf_file, chromosome, start_pos, K_in):

    # M_in, C_in, K_in = (MD_vars, MD_locs), (CIGAR_vars, CIGAR_locs), (variants, var_locs)
    # M, C, K = (MD_vars, MD_locs), (CIGAR_vars, CIGAR_locs), (variants, var_locs)
    K_i = 0
    K = (K_in[0], K_in[1] + (max(K_in[1]) +1000,))
    comb_err = ""
    new_vars_string = ""

    while K_i < len(K[0]):


        new_variant = K[0][K_i]
        K_Flag = "True"


        #new_vars_string += vcf_print(vcf_file, new_variant, K_Flag) # vcf_write(vcf_file, K[0][K_i], "True")
        new_vars_string += vcf_write(vcf_file, new_variant, K_Flag)
        K_i += 1

    return new_vars_string[1:], comb_err

def vcf_print(vcf_file, variant, TGFlag):

    variant_list = variant.split(",")
    vcf_line = variant_list[0].split(":") # chromosome and position
    vcf_line += [".", variant_list[1], variant_list[2]] # ID, ref, alt
    vcf_line += [".", "PASS"] # QUAL, FILTER

    new_vars_string = "chr" + variant_list[0] + "," + variant_list[1] + "," + variant_list[2]

    info = ["1000G="+TGFlag]

    vcf_line.append(";".join(info))

    print("\t".join(vcf_line))

    return ";"+new_vars_string


def vcf_write(vcf_file, variant, TGFlag):

    variant_list = variant.split(",")
    vcf_line = variant_list[0].split(":") # chromosome and position
    vcf_line += [".", variant_list[1].upper(), variant_list[2].upper()] # ID, ref, alt
    vcf_line += [".", "PASS"] # QUAL, FILTER

    new_vars_string = "chr" + variant_list[0] + "," + variant_list[1] + "," + variant_list[2]

    info = ["1000G="+TGFlag]

    vcf_line.append(";".join(info))

    #vcf_file.write("\t".join(vcf_line) + "\n")
    return ";"+new_vars_string

def write_CIGAR(Vs, Ls, S, length, CS):
    """
    This writes a new CIGAR string based on the known thosuand genomes variants and the existing CIGAR string

    Parameters
    ----------
    Vs : tuple
        The list of variants that exist in the read.
    Ls : tuple
        The list of locations of the variants in the read in the personal genome sequence.
    S : int
        The start coordinate for this read in the personal genome sequence.
    length : int
        The length of the read.
    CS : str
        The original CIGAR string.

    Returns
    -------
    C : str
        The new CIGAR string.

    """
    start_adjusted = S

    L_frag_i, L_frag_v = [], []
    adjustment = 0
    for i in range(len(Vs)):
        temp = Vs[i].split(",")
        if len(temp[1]) + len(temp[2]) != 2:
            L_frag_v.append(temp)
            if len(L_frag_i) > 0:
                L_frag_i.append(max(int(temp[0][temp[0].index(":")+1:])-S-adjustment, max(L_frag_i)))
            else:
                L_frag_i.append(int(temp[0][temp[0].index(":")+1:])-S-adjustment)
            adjustment += max(len(temp[1]) - len(temp[2]), 0)
    L_frag_i.append(length)


    if len(L_frag_v) == 0:
        return CS, start_adjusted

    else:

        CS_list = [[], []]
        val_temp = ""
        for i in CS:
            if i in ["M", "I", "D"]:
                CS_list[0].append(int(val_temp))
                CS_list[1].append(i)
                val_temp = ""
            else:
                val_temp += i
        #CS_list[0].append(max(L_frag_i))

        CS_k_list = [[], []]
        K_t = 0
        for i in range(len(L_frag_v)):
            if L_frag_i[i] >= 0:
                CS_k_list[0].append(L_frag_i[i] - K_t + 1)
                CS_k_list[1].append("M")
                K_t += L_frag_i[i] - K_t + 1
                if len(L_frag_v[i][2])>len(L_frag_v[i][1]):
                    CS_k_list[0].append(len(L_frag_v[i][2])-len(L_frag_v[i][1]))
                    CS_k_list[1].append("I")
                else:
                    CS_k_list[0].append(len(L_frag_v[i][1])-len(L_frag_v[i][2]))
                    CS_k_list[1].append("D")
            else:
                if len(L_frag_v[i][2])>len(L_frag_v[i][1]):
                    CS_k_list[0].append(len(L_frag_v[i][2])-len(L_frag_v[i][1]) + L_frag_i[i] + 1)
                    CS_k_list[1].append("I")
                else:
                    CS_k_list[0].append(len(L_frag_v[i][1])-len(L_frag_v[i][2]) + L_frag_i[i] + 1)
                    CS_k_list[1].append("D")

        C_i, K_i, C_tot, K_tot, total_dist = 0, 0, 0, 0, 0
        C_temp, K_temp = 0, 0
        CS_new = [[], []]

        while C_i < len(CS_list[0])-1 and K_i < len(CS_k_list[0]):

            # If the next K index is less than the next C index (paste K)
            if K_tot + CS_k_list[0][K_i] <= C_tot + CS_list[0][C_i]:
                CS_new[1].append(CS_k_list[1][K_i])
                CS_new[0].append(max(CS_k_list[0][K_i] - C_temp, 0))
                if CS_new[1][-1] != "D":
                    K_tot += CS_k_list[0][K_i]
                    K_temp += CS_k_list[0][K_i]
                    total_dist += CS_k_list[0][K_i]
                K_i += 1
                K_temp -= C_temp
                C_temp = 0
                if K_i < len(CS_k_list[1]) and CS_k_list[1][K_i] == "D":
                    CS_new[0].append(CS_k_list[0][K_i])
                    CS_new[1].append(CS_k_list[1][K_i])
                    K_i += 1
                if K_temp < 0:
                    C_temp -= K_temp
                    K_temp = 0
            elif C_i < len(CS_list[0]) and K_tot + CS_k_list[0][K_i] <= C_tot + CS_list[0][C_i] + CS_list[0][C_i+1] and CS_k_list[1][K_i] != "M":
                CS_new[1].append(CS_k_list[1][K_i])
                CS_new[0].append(max(CS_k_list[0][K_i] - C_temp, 0))
                if CS_new[1][-1] != "D":
                    K_tot += CS_k_list[0][K_i]
                    K_temp += CS_k_list[0][K_i]
                    total_dist += CS_k_list[0][K_i]
                K_i += 1
                K_temp -= C_temp
                C_temp = 0
                if K_i < len(CS_k_list[1]) and CS_k_list[1][K_i] == "D":
                    CS_new[0].append(CS_k_list[0][K_i])
                    CS_new[1].append(CS_k_list[1][K_i])
                    K_i += 1
                if K_temp < 0:
                    C_temp -= K_temp
                    K_temp = 0

            # If the next C index is less than the next K index (paste C)
            else:
                CS_new[1].append(CS_list[1][C_i])
                CS_new[0].append(max(CS_list[0][C_i] - K_temp, 0))
                if CS_new[1][-1] != "D":
                    C_tot += CS_list[0][C_i]
                    C_temp += CS_list[0][C_i]
                    total_dist += CS_k_list[0][K_i]
                C_i += 1
                C_temp -= K_temp
                K_temp = 0
                if C_i < len(CS_list[1]) and CS_list[1][C_i] == "D":
                    CS_new[0].append(CS_list[0][C_i])
                    CS_new[1].append(CS_list[1][C_i])
                    C_temp += CS_list[0][C_i]
                    C_i += 1




        while K_i < len(CS_k_list[0]):
            CS_new[0].append(max(CS_k_list[0][K_i] - C_temp, 0))
            CS_new[1].append(CS_k_list[1][K_i])
            if CS_new[1][-1] != "D":
                K_tot += CS_k_list[0][K_i]
                K_temp += CS_k_list[0][K_i]
            K_i += 1
            K_temp -= C_temp
            C_temp = 0

        while C_i < len(CS_list[0]):
            CS_new[0].append(max(CS_list[0][C_i] - K_temp, 0))
            CS_new[1].append(CS_list[1][C_i])
            if CS_new[1][-1] != "D":
                C_tot += CS_list[0][C_i]
                C_temp += CS_list[0][C_i]
            C_i += 1
            C_temp -= K_temp
            K_temp = 0


        n_vals = []
        n_labs = []
        tot = 0
        for i in range(len(CS_new[0])):
            if CS_new[0][i] != 0:
                if len(n_labs) > 0 and CS_new[0][i] == n_labs[-1]:
                    n_vals[-1] += CS_new[0][i]
                    tot += n_vals[-1]
                elif tot + CS_new[0][i] < length or CS_new[1][i] == "D":
                    n_vals.append(CS_new[0][i])
                    n_labs.append(CS_new[1][i])
                    if n_labs[-1] != "D": tot += n_vals[-1]
                else:
                    n_vals.append(min(length-tot, length))
                    n_labs.append(CS_new[1][i])
                    if n_labs[-1] != "D": tot = length



        n_vals2 = [n_vals[0]]
        n_labs2 = [n_labs[0]]
        for i in range(len(n_vals)-1):
            if n_labs[i+1] == n_labs2[-1]:
                n_vals2[-1] += n_vals[i+1]
            else:
                n_labs2.append(n_labs[i+1])
                n_vals2.append(n_vals[i+1])

        C = ""
        for i in range(len(n_vals2)):
            C += str(n_vals2[i])
            C += n_labs2[i]

        return C, start_adjusted

# %% Scan SAM file debug

sam_file = open(args.output, 'w')
vcf_file = 0 #vcf_file = open(vcf_name, 'w')
#errors_file = open(errors_name, 'w')
sam_file.write(file_headers.sam_bowtie_hg38())
sam_file.close()
sam_file = open(args.output, 'a')

for line_og in open(args.input):
 line = line_og[:-1].split("\t")
 if line[0][0] != "@":
  flag = str(bin(int(line[1])))
  # Check that it mapped to a sequence
  if (line[2] == "chrM" or line[2] == "chrY") and line[6] == "=":

    #sam_file.write("\t".join(line) + "\n")
    pass


  elif line[2] != "*" and line[6] == "=" and line[1] in ["99", "147", "83", "163"]: # flag[-6:-4] in ["10", "01"]:

    # Create the new output line
    new_sam = line.copy()

    # Pull the important values from the line
    sequence = line[9] # sequence of the read
    seq_len = len(sequence) # length of the sequence
    frag_len = abs(int(line[8])) # length of the fragment
    start_ind = int(line[3]) # start index within the personal genome sequence
    mate_start_ind = int(line[7]) # start index of the mate within the personal genome seuence

    # Re-organize the known variants
    seq_name, variants_short, var_locs_short, variants, var_locs = find_name_and_vars(line[2], min(start_ind, mate_start_ind), frag_len)
    variants_long, var_locs_long = variants, var_locs
    chromosome = seq_name[:seq_name.index(":")]
    start_pos = int(seq_name[seq_name.index(":")+1:seq_name.index("-")])

    seq_name, variants_short, var_locs_short, variants, var_locs = find_name_and_vars(line[2], start_ind, seq_len)

    new_vars_string = ""
    if len(variants_short) > 0:
        new_vars_string, comb_err = process_variants(vcf_file, chromosome, start_pos,
                                                     (variants_short, var_locs_short))

    if len(new_vars_string) > 0:
        new_sam.append("rv:Z:"+new_vars_string)

    seq_name, variants_short, var_locs_short, variants, var_locs = find_name_and_vars(line[2], min(start_ind, mate_start_ind), frag_len)

    new_vars_string = ""
    if len(variants_short) > 0:
        new_vars_string, comb_err = process_variants(vcf_file, chromosome, start_pos,
                                                     (variants_short, var_locs_short))

    if len(new_vars_string) > 0:
        new_sam.append("fv:Z:"+new_vars_string)

    # Find the start coordinate of this read and its mate
    start_coord = int(seq_name[seq_name.index(":")+1:seq_name.index("-")]) + int(line[3])
    mate_coord = int(seq_name[seq_name.index(":")+1:seq_name.index("-")]) + int(line[7])
    new_CIGAR = ""
    for i in range(len(var_locs_long)):
        var = variants_long[i].split(",")
        if var_locs_long[i] < int(line[3]):
            start_coord = max(int(var[0][var[0].index(":")+1:]) + int(line[3]) - var_locs_long[i] - len(var[2]) + len(var[1]) + 1, int(var[0][var[0].index(":")+1:])+1)
            if len(var[1]) != 1 and int(var[0][var[0].index(":")+1:])+len(var[1])-1 > start_coord:
                # This is a deletion
                new_CIGAR = str(start_coord-(int(var[0][var[0].index(":")+1:])+len(var[1])-1)) + "D"
            if len(var[2]) != 1 and int(var[0][var[0].index(":")+1:])+len(var[2])-1 > start_coord:
                # This is an insertion
                new_CIGAR = str(start_coord-(int(var[0][var[0].index(":")+1:])+len(var[2])-1)) + "I"
        if var_locs_long[i] < int(line[7]):
            mate_coord = int(var[0][var[0].index(":")+1:]) + int(line[7]) - var_locs_long[i] - len(var[2]) + len(var[1]) + 1

    # Create a new CIGAR string
    seq_name, variants_short, var_locs_short, variants, var_locs = find_name_and_vars(line[2], start_ind, seq_len)
    new_CIGAR, start_coord = write_CIGAR(variants_short, var_locs_short, start_coord, seq_len, line[5])

    new_sam[2] = chromosome # chromosome
    new_sam[3] = str(start_coord) # position of read
    new_sam[5] = new_CIGAR # CIGAR string
    new_sam[7] = str(mate_coord) # position of mate
    if line[1] in ["99", "163"] and "-" in new_sam[8]:
        new_sam[8] = new_sam[8][1:]
    sam_file.write("\t".join(new_sam) + "\n")

sam_file.close()
