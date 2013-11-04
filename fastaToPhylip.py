#! /usr/bin/python

# Script to convert an alignment in the fasta format to the phylip format.
# The output format is sequential.
#
# Usage:
# ------
#
# fastaToPhylip.py inputFile



# import
# ------

import sys

from Bio import SeqIO



# main
# ----



# parse the arguments

input_file = sys.argv[1]

output_file = input_file.split(".fasta")[0] + ".phylip"



# open the input file

fi = open(input_file, "r")

input_iter = SeqIO.parse(fi, format = "fasta")



# open the output file

fo = open(output_file, "w")



# go through the input file once to get the sequence names and lengths, and the
# total number of sequences

n = 0

seq_names = []

seq_lengths = []

for seq in input_iter :
    
    n += 1
    
    seq_names.append(seq.name)
    
    seq_lengths.append(len(str(seq.seq)))
    
# check that all sequences are the same length

assert len(set(seq_lengths)) == 1
    
# get the length of the longest name

len_names = [len(x) for x in seq_names]

len_names_max = max(len_names)



# reopen the input file

fi.close()

fi = open(input_file, "r")

input_iter = SeqIO.parse(fi, format = "fasta")



# go through the sequences and write them to the output file

fo.write(str(n) + " " + str(seq_lengths[0]) + "\n")

for seq in input_iter :
    
    fo.write(seq.name + " " * (len_names_max - len(seq.name)) + "  " + str(seq.seq) + "\n")

fo.close()

fi.close()
    

    
    
