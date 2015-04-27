### * Description

# Script to get unique fasta sequences from a list of input sequences.
#
# Output is the unique sequences in a fasta format, using a hash of the
# sequence as the sequence name.
#
# In addition to the Python documentation, useful information can be found here
# for how to use stdin as the input stream:
# http://stackoverflow.com/questions/7576525/optional-stdin-in-python-with-argparse
# For how to add the defaults in help message:
# http://stackoverflow.com/questions/12151306/argparse-way-to-include-default-values-in-help

### ** Usage

# python script.py -h

### * Setup

### ** Import
import sys
import hashlib
import argparse
from Bio import SeqIO

### ** Argument parser
parser = argparse.ArgumentParser(
    description =
    "Report unique fasta sequences from the input. Output is in fasta format, \n"
    "with sequence names being a hash of the actual sequence.",
    epilog =
    "Example usage:\n"
    "--------------\n"
    "python " + sys.argv[0] + " mySeq.fasta -o myUniqSeq.fasta\n"
    "cat *.fasta | python " + sys.argv[0] + " > myUniqSeq.fasta",
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input", metavar = "INPUT",
                    nargs = "?",
                    help = "A fasta file. If no file is provided, stdin is used",
                    type = argparse.FileType("r"),
                    default = sys.stdin)
parser.add_argument("-o", "--output", 
                    help = "Output file. If no file is provided, stdout is used",
                    default = False)

### * Functions

### ** hashSequence(sequence)

def hashSequence(sequence) :
    """
    TAKES
    sequence: a string
    RETURNS
    The sha512 hash for it
    """
    m = hashlib.sha512()
    m.update(sequence)
    return(m.hexdigest())

### * Run

if (__name__ == "__main__") :
    args = parser.parse_args()

### ** Open output
if not args.output :
    fo = sys.stdout
else :
    fo = open(args.output, "w")
    
### ** Parse the input
o = set()
for seq in SeqIO.parse(args.input, "fasta") :
    sequence = str(seq.seq)
    hashValue = hashSequence(sequence)
    if (hashValue not in o) :
        fo.write(">" + hashValue + "\n")
        fo.write(sequence + "\n")
        o.add(hashValue)

### * End
sys.stdin.close()
fo.close()
