#! /usr/bin/python

# restrictionCut.py

# -----------
# description
# -----------

# This script performs restriction cuts on a fasta file and produce a fasta file 
# containing the sequences and some graphs.

# NOTE: the cutSequence function proceeds by reading the sequence and cutting
# the next occurrence of the restriction site. In a case such as a "CCC-CC"
# restriction site and a "ACCCCCCG" sequence, which could produce in reality
# both "ACCC CCCG" and "ACCCC CCG" cuts, only the first one will be produced
# by the function.




# -------
# imports
# -------



import sys

import textwrap

import pickle

from Bio import SeqIO

from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord



# -------
# classes
# -------



class openFastaFile :


    
    """
    This class helps to handle fasta files.
    
    """



    def __init__(self, filename) :
        
        """
        Opens a fasta file
        
        Takes
        filename: the path to the file
        
        """
        
        # open the file
        
        self.lastOpenedFile = filename
        
        self.fileHandle = open(filename, 'r')
        
        self.fastaFileHandle = SeqIO.parse(self.fileHandle, 'fasta')


    
    def __iter__(self) :
        
        """
        Iterator to go through the fasta file.
        
        """
        
        return self



    def reset(self) :

        """
        Closes and opens the file.
        
        """
        
        self.close()
        
        self.__init__(self.lastOpenedFile)



    def next(self) :

        """
        Returns the next element from the file.
        
        """
        
        return self.fastaFileHandle.next()



    def close(self) :
        
        """
        Closes the opened file.
        
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
    
    """
    Takes a Biopython sequence and a motif and return the sequences cut at 
    this motif.
    
    Takes
    sequence: the sequence as a string, a Seq or a SeqRecord
    motif: a string containing a '-' at the cut site; e.g. 'g-aattc' for EcoRI
    
    Returns
    a list of Biopython sequences SeqRecord, with names corresponding to the 
    input sequence, the cut and the location in the input sequence
    
    # NOTE: the cutSequence function proceeds by reading the sequence and cutting
    # the next occurrence of the restriction site. In a case such as a "CCC-CC"
    # restriction site and a "ACCCCCCG" sequence, which could produce in reality
    # both "ACCC CCCG" and "ACCCC CCG" cuts, only the first one will be produced
    # by the function.
    
    """
    
    
    # clean the cutting motif
    
    pattern = motif.replace('-', '').upper()
    
    cutSite = motif.find('-')
    
    if (cutSite == -1) :
    
        cutSite = len(motif) # if no '-' is found, the cut is at the end of the motif
    
    # extract the sequence string
    
    try :
    
        sequenceString = sequence.seq.upper() # if sequence is a SeqRecord
    
        name = sequence.name+'_'+motif+'_'
    
    except :
    
        try :
    
            sequenceString = sequence.tostring().upper() # if sequence is a Seq
    
            name = motif+'_'
    
        except :
    
            sequenceString = sequence.upper() # sequence should be a string
    
            name = motif+'_'
    
    # find the first occurrence
    
    fragments = []
    
    nextSite = sequenceString.find(pattern)
    
    lastCut = 1
    
    lastPosition = len(sequenceString)
    
    # loop
    
    while (nextSite>0) :
    
        # while a site is found
    
        newFragment = SeqRecord(sequenceString[ : nextSite+cutSite])
    
        newFragment.name = (name + str(lastCut) + '_' + 
                           str(lastCut + nextSite + cutSite - 1))
    
        lastCut = lastCut + nextSite + cutSite
    
        sequenceString = sequenceString[nextSite + cutSite : ]
        
        fragments.append(newFragment)
        
        nextSite = sequenceString.find(pattern)
        
    # add the remaining sequence
    
    if (sequenceString != '') :
    
        name = name + str(lastCut) + '_' + str(lastPosition)
        
        lastFragment = SeqRecord(sequenceString)
        
        lastFragment.name = name
        
        fragments.append(lastFragment)
        
    # return
    
    return fragments



def saveSequenceList(sequenceList, filename, mode='w') :
    
    """
    Saves a list of Biopython sequences SeqRecord to a fasta file.
    
    Takes
    sequenceList: list of SeqRecord
    filename: the name of the file
    mode: the mode of file opening
    
    Returns
    none
    
    """
    
    f=open(filename, mode)
    
    for sequence in sequenceList :
        
        f.write('>'+sequence.name+'\n')
        
        f.write(sequence.seq + "\n")
        
    f.close()



def getCutFragmentLengths(sequence, motif) :
    """
    Takes a sequence and a motif and returns a list of the length of the 
    produced fragments.
    
    Takes
    sequence: the sequence
    motif: a string containing a '-' at the cut site ; e.g. 'g-aattc' for EcoRI
    
    Returns
    a list of fragment lengths
    
    """
    
    fragments = cutSequence(sequence, motif)
    
    lengths = [len(x) for x in fragments]
    
    return lengths
    
    

def cutSequenceList(sequenceList, motif) :
    
    """
    cutSequence() for a list. Concatenates the results into a list.
    
    """
    
    fragments = []
    
    for sequence in sequenceList :
        
        newFragments = cutSequence(sequence, motif)
        
        fragments += newFragments
        
    return fragments
    


def getCutFragmentLengthsList(sequenceList, motif) :
    
    """
    getCutFragmentLengths() for a lsit. Concatenates all the lengths into a list.
    
    """
    
    lengths = []
    
    for sequence in sequenceList :
        
        newLengths = cutFragmentLengths(sequence, motif)
        
        lengths += newLengths
        
    return lengths



def cutSequenceFile(filename, motif, outputFilename = '') :
    
    """
    cutSequence() on the sequences contained in a fasta file. If an output 
    filename is given, saves all the results to the file. If not, returns a 
    list.
    
    NOTE: it appears that using an outputFilename makes it slow. It is
    better to do it all in memory, store the results as a list and then
    only write it to the disk.
    
    """
    
    f = open(filename, 'r')
    
    inputFile = SeqIO.parse(f, 'fasta')
    
    if (outputFilename == '') :
        
        # if no filename is given
        
        fragments = []
        
        for sequence in inputFile :
            
            print("Cutting " + sequence.name + ", length " + 
                  str(len(sequence)) + " bp")
            
            newFragments = cutSequence(sequence, motif)
            
            print(str(len(newFragments)) + " fragments generated")
            
            fragments += newFragments
            
        f.close()
        
        return fragments
        
    else :
        
        # if a filename is given
        
        f = open(outputFilename, "w")
        
        f.close() # just to erase its content
    
        for sequence in inputFile :
            
            print("Cutting " + sequence.name + ", length " + 
                  str(len(sequence)) + " bp")
            
            fragments=cutSequence(sequence, motif)
            
            print(str(len(fragments)) + " fragments generated")
            
            saveSequenceList(fragments, outputFilename, 'a')
            
        f.close()



def getCutFragmentLengthsFile(filename, motif) :
    
    """
    getCutFragmentLengths() on the sequences from a fasta file. Returns a list.
    
    """
    
    lengths = []
    
    f = open(filename, 'r')
    
    inputFile = SeqIO.parse(f, 'fasta')
    
    for sequence in inputFile :
    
        print("Cutting " + sequence.name)
    
        lengths += cutFragmentLengths(sequence, motif)
        
    f.close()
    
    return lengths



def filterFastaFile(inputFilename, outputFilename, patternToFind) :
    
    """
    Takes only the input sequences containing a given pattern in their name.
    Writes them to a fasta file.
    
    (case sensitive)
    
    Takes
    inputFilename
    outputFilename
    patternToFind
    
    Returns
    none
    
    """
    
    inputFile = open(inputFilename, 'r')
    
    outputFile = open(outputFilename, 'w')
    
    sequences = SeqIO.parse(inputFile, 'fasta')
    
    for sequence in sequences :
    
        if (patternToFind in sequence.name) :
            
            outputFile.write('>' + sequence.name + '\n')
            
            outputFile.write(sequence.seq.data + '\n')
            
    inputFile.close()
    
    outputFile.close()



def saveListToFile(inputList, outputFilename) :

    """
    Writes a numerical list to a file.
    
    Takes
    inputList
    outputFilename
    
    Returns
    none
    
    """
    
    f = open(outputFilename, 'w')
    
    for element in inputList :
    
        f.write(str(element) + '\n')
    
    f.close()



def sizeFilter(sequenceList, minSize = 1, maxSize = 0) :
    
    """
    Filter a list of sequences and only keep the ones comprised between
    the min and max sizes, included.
    
    Takes
    minSize
    maxSize: 0 if no max size
    
    Returns
    a list of sequences
    
    """
    
    output = []
    
    if (maxSize > 0) :
    
        for seq in sequenceList :
        
            if (len(seq) >= minSize and len(seq) <= maxSize) :
                
                output.append(seq)
                
    else :
        
        for seq in sequenceList :
            
            if (len(seq) >= minSize) :
                
                output.append(seq)
                
    return output
    
        

#~ # ------
#~ # script
#~ # ------
#~ 
#~ if (__name__=='__main__') :
#~ 
    #~ # input file
    #~ inputFileFull='Gasterosteus_aculeatus.BROADS1.56.dna.toplevel.fa'
    #~ inputFile2='Ga_chromosomesOnly.fa'
    #~ inputFile=inputFile2
#~ 
    #~ # motives
    #~ motives=['ctgca-g', 'g-aattc', 'cctgca-gg']
#~ 
    #~ # do the cuts
    #~ for motif in motives :
        #~ lengths=cutFragmentLengthsFile(inputFile, motif)
        #~ # save the lengths
        #~ f=open(inputFile+'.'+motif, 'w')
        #~ pickle.dump(lengths, f)
        #~ f.close()
