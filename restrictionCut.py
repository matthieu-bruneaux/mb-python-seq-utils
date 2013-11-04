#! /usr/bin/python

# restrictionCut.py

# -----------
# description
# -----------

# This script performs restriction cuts on a fasta file and produce a fasta file containing the sequences and some graphs.


# -------
# imports
# -------

import sys
import pickle
import toolBox
from Bio import SeqIO
from Bio.Seq import Seq, SeqRecord


# -------
# classes
# -------

class openFastaFile :
    """This class helps to hanlde fasta files.
    """

    def __init__(self, filename) :
        """Opens a fasta file
        Input:
        filename = the path to the file
        """
        # open the file
        self.lastOpenedFile=filename
        self.fileHandle=open(filename, 'r')
        self.fastaFileHandle=SeqIO.parse(self.fileHandle, 'fasta')

    
    def __iter__(self) :
        """Iterator for going through the fasta file.
        """
        return self.fastaFileHandle

    def reset(self) :
        """Closes and opens the file.
        """
        self.close()
        self.__init__(self.lastOpenedFile)

    def next(self) :
        """Returns the next element from the file.
        """
        return self.fastaFileHandle.next()

    def close(self) :
        """Closes the opened file.
        """
        self.fileHandle.close()
        try :
            del(self.fastaFileHandle)
        except :
            pass


# ---------
# functions
# ---------

def cutSequence(sequence, motif) :
    """Takes a Biopython sequence and a motif and return the sequences cut at this motif.
    Input:
    sequence = the sequence
    motif = a string containing a '-' at the cut site ; e.g. 'g-aattc' for EcoRI
    Output:
    a list of Biopython sequences SeqRecord, with names corresponding to the input sequence, the cut and the location in the input sequence
    """
    # clean the cutting motif
    pattern=motif.replace('-', '').upper()
    cutSite=motif.find('-')
    if (cutSite==-1) :
        cutSite=len(motif) # if no '-' is found the cut is at the end of the motif
    # extract the sequence string
    try :
        sequenceString=sequence.seq.data.upper() # if sequence is a SeqRecord
        name=sequence.name+'_'+motif+'_'
    except :
        try :
            sequenceString=sequence.data.upper() # if sequence is a Seq
            name=motif+'_'
        except :
            sequenceString=sequence.upper() # sequence should be a string
            name=motif+'_'
    # find the first occurrence
    fragments=[]
    nextSite=sequenceString.find(pattern)
    lastCut=1
    lastPosition=len(sequenceString)
    # loop
    while (nextSite>0) :
        # while a site is found
        newFragment=SeqRecord(sequenceString[:nextSite+cutSite])
        newFragment.name=name+str(lastCut)+'_'+str(lastCut+nextSite+cutSite-1)
        lastCut=lastCut+nextSite+cutSite
        sequenceString=sequenceString[nextSite+cutSite:]
        fragments.append(newFragment)
        nextSite=sequenceString.find(pattern)
    # add the remaining sequence
    if (sequenceString!='') :
        name=name+str(lastCut)+'_'+str(lastPosition)
        lastFragment=SeqRecord(sequenceString)
        lastFragment.name=name
        fragments.append(lastFragment)
    # return
    return fragments


def saveSequenceList(sequenceList, filename, mode='w') :
    """Saves a list of Biopython sequences SeqRecord to a fasta file.
    Input:
    sequenceList = list
    filename = the name of the file
    mode = the mode of file opening
    Output:
    none
    """
    f=open(filename, mode)
    for sequence in sequenceList :
        f.write('>'+sequence.name+'\n')
        f.write(sequence.seq+'\n')
    f.close()


def cutFragmentLengths(sequence, motif) :
    """Takes a sequence and a motif and returns a list of the length of the produced fragments.
    Input:
    sequence = the sequence
    motif = a string containing a '-' at the cut site ; e.g. 'g-aattc' for EcoRI
    Output:
    a list of fragment lengths
    """
    fragments=cutSequence(sequence, motif)
    lengths=[len(x) for x in fragments]
    return lengths


def cutSequenceList(sequenceList, motif) :
    """cutSequence() for a list. Concatenates the results in a list.
    """
    fragments=[]
    for sequence in sequenceList :
        newFragments=cutSequence(sequence, motif)
        fragments+=newFragments
    return fragments


def cutFragmentLengthsList(sequenceList, motif) :
    """cutFragmentLengths() for a lsit. Concatenates all the lengths in a list.
    """
    lengths=[]
    for sequence in sequenceList :
        newLengths=cutFragmentLengths(sequence, motif)
        lengths+=newLengths
    return lengths


def cutSequenceFile(filename, motif, outputFilename='') :
    """cutSequence() on the sequences contained in a fasta file. If an output filename is given, saves all the results to the file. If not, returns a list.
    """
    f=open(filename, 'r')
    inputFile=SeqIO.parse(f, 'fasta')
    if (filename=='') :
        # if no filename is given
        fragments=[]
        for sequence in inputFile :
            print sequence.name
            newFragments=cutSequence(sequence, motif)
            fragments+=newFragments
        f.close()
        return fragments
    else :
        # if a filename is given
        for sequence in inputFile :
            print sequence.name
            fragments=cutSequence(sequence, motif)
            saveSequenceList(fragments, outputFilename, 'a')
        f.close()


def cutFragmentLengthsFile(filename, motif) :
    """cutFragmentLengths() on the sequences from a fasta file. Returns a list.
    """
    lengths=[]
    f=open(filename, 'r')
    inputFile=SeqIO.parse(f, 'fasta')
    for sequence in inputFile :
        print sequence.name
        lengths+=cutFragmentLengths(sequence, motif)
    f.close()
    return lengths


def filterFastaFile(inputFilename, outputFilename, patternToFind) :
    """Takes only the input sequences containing a given pattern in their name.
    Writes them to a fasta file.
    Input:
    inputFilename
    outputFilename
    patternToFind
    Output:
    none
    """
    inputFile=open(inputFilename, 'r')
    outputFile=open(outputFilename, 'w')
    sequences=SeqIO.parse(inputFile, 'fasta')
    for sequence in sequences :
        if (patternToFind in sequence.name) :
            outputFile.write('>'+sequence.name+'\n')
            outputFile.write(sequence.seq.data+'\n')
    inputFile.close()
    outputFile.close()


def saveListToFile(inputList, outputFilename) :
    """Writes a numerical list to a file.
    Input:
    inputList
    outputFilename
    Output:
    none
    """
    f=open(outputFilename, 'w')
    for element in inputList :
        f.write(str(element)+'\n')
    f.close()


# ------
# script
# ------

if (__name__=='__main__') :

    # input file
    inputFileFull='Gasterosteus_aculeatus.BROADS1.56.dna.toplevel.fa'
    inputFile2='Ga_chromosomesOnly.fa'
    inputFile=inputFile2

    # motives
    motives=['ctgca-g', 'g-aattc', 'cctgca-gg']

    # do the cuts
    for motif in motives :
        lengths=cutFragmentLengthsFile(inputFile, motif)
        # save the lengths
        f=open(inputFile+'.'+motif+'.'+toolBox.timeStamp(), 'w')
        pickle.dump(lengths, f)
        f.close()
