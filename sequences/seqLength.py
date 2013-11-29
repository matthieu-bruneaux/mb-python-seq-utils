#! /usr/bin/python

# script to produce a list of the sequences length from a fasta file

import sys
from Bio import SeqIO

# parse the input
inputFile=sys.argv[1]
outputFile=sys.argv[1]+'.seqLength'

# go through the file
f=open(inputFile, 'r')
f2=open(outputFile, 'w')
sequences=SeqIO.parse(f, 'fasta')

for seq in sequences :
    f2.write(seq.name+'\t'+str(len(seq))+'\n')
    print(seq.name+'\t'+str(len(seq))+'\n')

f2.close()
f.close()
