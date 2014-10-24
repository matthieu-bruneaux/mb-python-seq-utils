# script to convert fasta files with line wrapping to fasta file with
# each sequence on one line

import sys


if len(sys.argv) < 3:
    print("> Enter input and output files as arguments.")
    sys.exit(1)

nameIn = sys.argv[1]
nameOut = sys.argv[2]


fIn=open(nameIn, 'r')
fOut=open(nameOut ,'w')

line1=fIn.readline()

line2=fIn.readline()

while (line2) :
    if (('>' in line2) | ('>' in line1)) :
        fOut.write(line1)
        print line1
    else :
        fOut.write(line1.rstrip())
    line1=line2
    line2=fIn.readline()

fOut.write(line1)

fIn.close()
fOut.close()
