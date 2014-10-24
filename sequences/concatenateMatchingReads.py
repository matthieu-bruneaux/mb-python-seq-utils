#! /usr/bin/python

# ------------------------------------------------------------------- #
# Concatenate matching reads in two input files (forward and reverse) #
# ------------------------------------------------------------------- #



# Format :
# - linux  : ./concatenateMatchingReads.py inputFileForward inputFileReverse outputFile
# - windows : python concatenateMatchingReads.py inputFileForward inputFileReverse outputFile

# Note :
# This script assumes that each file has exactly 4 lines per sequences, with no blank lines at all
# (not at the end either), so that each (4*n + 1) is a new sequence.


# import
# ------



import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



# define a function from Biopython not defined in my old Biopython
# ----------------------------------------------------------------



def doReverseComplement(inputSeqRecord, id=False, name=False, description=False, features=True, annotations=False, letter_annotations=True, dbxrefs=False) :

    from Bio.SeqRecord import SeqRecord    
    from Bio.Seq import MutableSeq #Lazy to avoid circular imports 
    if isinstance(inputSeqRecord.seq, MutableSeq): 
        #Currently the MutableSeq reverse complement is in situ 
        answer = SeqRecord(inputSeqRecord.seq.toseq().reverse_complement()) 
    else: 
        answer = SeqRecord(inputSeqRecord.seq.reverse_complement()) 
    if isinstance(id, basestring): 
        answer.id = id 
    elif id: 
        answer.id = inputSeqRecord.id 
    if isinstance(name, basestring): 
        answer.name = name 
    elif name: 
        answer.name = inputSeqRecord.name 
    if isinstance(description, basestring): 
        answer.description = description 
    elif description: 
        answer.description = inputSeqRecord.description 
    if isinstance(dbxrefs, list): 
        answer.dbxrefs = dbxrefs 
    elif dbxrefs: 
        #Copy the old dbxrefs 
        answer.dbxrefs = inputSeqRecord.dbxrefs[:] 
    if isinstance(features, list): 
        answer.features = features 
    elif features: 
        #Copy the old features, adjusting location and string 
        l = len(answer) 
        answer.features = [f._flip(l) for f in inputSeqRecord.features] 
        #The old list should have been sorted by start location, 
        #reversing it will leave it sorted by what is now the end position, 
        #so we need to resort in case of overlapping features. 
        #NOTE - In the common case of gene before CDS (and similar) with 
        #the exact same locations, this will still maintain gene before CDS 
        answer.features.sort(key=lambda x : x.location.start.position) 
    if isinstance(annotations, dict): 
        answer.annotations = annotations 
    elif annotations: 
        #Copy the old annotations, 
        answer.annotations = inputSeqRecord.annotations.copy() 
    if isinstance(letter_annotations, dict): 
        answer.letter_annotations = letter_annotations 
    elif letter_annotations: 
        #Copy the old per letter annotations, reversing them 
        for key, value in inputSeqRecord.letter_annotations.iteritems(): 
            answer._per_letter_annotations[key] = value[::-1] 
    return answer 



# define a function to concatenate two sequences with fastq quality information
# -----------------------------------------------------------------------------



def concatenateSequences(seq1, seq2, name) :

    # TAKES
    # seq1, seq2 : the sequences to be concatenated (SeqRecord)
    # name : the name of the output SeqRecord

    cId = name
        
    cSeq = Seq(seq1.seq.data + seq2.seq.data)
    
    cPhredQuality = seq1.letter_annotations["phred_quality"] + seq2.letter_annotations["phred_quality"]
    
    cLAnnotations = {"phred_quality" : cPhredQuality}
    
    concatenatedSequence = SeqRecord(cSeq, id = cId, letter_annotations = cLAnnotations, description = "")
    
    return concatenatedSequence



# parameters
# ----------



inputFileFormat = "fastq-sanger"
outputFileFormat = "fastq-sanger"

inputForwardTempFile = "input.forward.temp.file"
inputReverseTempFile = "input.reverse.temp.file"

counterStep = 100

# parse the input
# ---------------



if (len(sys.argv) >= 3) :
    
    inputFileNameForward = sys.argv[1]
    inputFileNameReverse = sys.argv[2]
    outputFileName = sys.argv[3]

else : # testing purpose

    inputFileNameForward = "./input/s_1_2_sequence_CACTCC_For_test.fastq"
    inputFileNameReverse = "./input/s_1_2_sequence_CACTCC_Rev_test.fastq"
    outputFileName = "./output"



# open the files
# --------------



inputForwardFileHandler = open(inputFileNameForward, "r")
inputReverseFileHandler = open(inputFileNameReverse, "r")
outputFileHandler = open(outputFileName, "w")



# go once through each input file and build a dictionary for each with
# the name of the sequence as key and the position in the file as a value
# -----------------------------------------------------------------------



# first define a handy function to read 4 lines at a time

def readNextSequence(fileHandler) :
    
    # TAKES
    # fileHandler : an open file handler
    
    # DOES
    # moves forward 4 lines in the file
    
    # RETURNS
    # a tuple containing the position of the file pointer before reading the lines
    # and the first line read
    # by moving with fileHandler.seek to this position, the same 4 lines are read
    # returns False when the end of the file is reached
    
    position = fileHandler.tell()
    
    firstLine = fileHandler.readline()
    
    for i in range(0, 3) :
    
        fileHandler.readline()
    
    if (firstLine == '') :
        
        return False
    
    else :
    
        return (position, firstLine)

    
    
# now read the input files

inputForwardDictionary = {}
inputReverseDictionary = {}



# forward input file

print "reading the sequence names from the forward input file..."

nextSequenceInformation = readNextSequence(inputForwardFileHandler)

while (nextSequenceInformation) :
    
    sequenceName = nextSequenceInformation[1][1 : -2]
    
    inputForwardDictionary[sequenceName] = nextSequenceInformation[0]

    nextSequenceInformation = readNextSequence(inputForwardFileHandler)


# reverse input file

print "reading the sequence names from the reverse input file..."


nextSequenceInformation = readNextSequence(inputReverseFileHandler)

while (nextSequenceInformation) :
    
    sequenceName = nextSequenceInformation[1][1 : ].rstrip("12\r\n")
    
    inputReverseDictionary[sequenceName] = nextSequenceInformation[0]

    nextSequenceInformation = readNextSequence(inputReverseFileHandler)



# transform the two lists of names to sets and determine their intersection
# prepare a list of names to be used to write the two new input files with
# their sequences ordered
# -------------------------------------------------------------------------



inputForwardNameSet = set(inputForwardDictionary.keys())
inputReverseNameSet = set(inputReverseDictionary.keys())

commonSequenceNameSet = inputForwardNameSet & inputReverseNameSet

commonSequenceNameList = list(commonSequenceNameSet)

print "forward input file : " + str(len(inputForwardNameSet)) + " sequences"
print "reverse input file : " + str(len(inputReverseNameSet)) + " sequences"
print "number of common names : " + str(len(commonSequenceNameSet)) +"\n"



# for each input file, use the ordered list of names and the dictionary to
# produce the temporary input file with the ordered common sequences
# ------------------------------------------------------------------------



# forward temp file

inputForwardTempFileHandler = open(inputForwardTempFile, "w")

for sequenceName in commonSequenceNameList :
    
    sequenceStart = inputForwardDictionary[sequenceName]
    
    inputForwardFileHandler.seek(sequenceStart, 0)
    
    for i in range(0, 4) :
        
        inputForwardTempFileHandler.write(inputForwardFileHandler.readline())
        
inputForwardTempFileHandler.close()


# reverse temp file

inputReverseTempFileHandler = open(inputReverseTempFile, "w")

for sequenceName in commonSequenceNameList :
    
    sequenceStart = inputReverseDictionary[sequenceName]
    
    inputReverseFileHandler.seek(sequenceStart, 0)
    
    for i in range(0, 4) :
        
        inputReverseTempFileHandler.write(inputReverseFileHandler.readline())
        
inputReverseTempFileHandler.close()



# go sequence by sequence through the two temporary input file, test if the
# names match as expected and write the final output
# -------------------------------------------------------------------------



print "Writing the concatenated common sequences... wait a little"

noSequence = 0

inputForwardTempFileHandler = open(inputForwardTempFile, "r")
inputReverseTempFileHandler = open(inputReverseTempFile, "r")

inputForwardTempFileIterator = SeqIO.parse(inputForwardTempFileHandler, format = inputFileFormat)
inputReverseTempFileIterator = SeqIO.parse(inputReverseTempFileHandler, format = inputFileFormat)


for (forwardCommonSequence, reverseCommonSequence) in zip(inputForwardTempFileIterator, inputReverseTempFileIterator) :
    
    # do the reverse complement
        
    reverseSeqRC = doReverseComplement(reverseCommonSequence)

    # build the concatenated sequence
    
    newName = forwardCommonSequence.name[0 : -1] + "1/2"
    
    concatenatedSequence = concatenateSequences(forwardCommonSequence, reverseSeqRC, newName)

    # write the sequence to the output file
    
    SeqIO.write([concatenatedSequence], outputFileHandler, outputFileFormat)

    # counter
    
    noSequence += 1
    
    if (noSequence % counterStep == 0) :
        
        print str(noSequence) + " sequences done"


print "\nforward input file : " + str(len(inputForwardNameSet)) + " sequences"
print "reverse input file : " + str(len(inputReverseNameSet)) + " sequences"
print "number of common names : " + str(len(commonSequenceNameSet)) +" - done -\n"



# close the files
# ---------------

inputForwardTempFileHandler.close()
inputReverseTempFileHandler.close()

os.remove(inputForwardTempFile)
os.remove(inputReverseTempFile)

inputForwardFileHandler.close()
inputReverseFileHandler.close()
outputFileHandler.close()
