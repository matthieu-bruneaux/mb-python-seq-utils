#! /usr/bin/python

# Extract the GO terms corresponding to UniProt identifiers.
# 
# Take two input files (list of UniProt identifiers and gene ontology 
# association file) and extract for each UniProt identifier all the 
# corresponding GO terms.
# The output file is a 4-column table with:
# * UniProt identifier
# * GO-BP
# * GO-CC
# * GO-MF
# The columns for the GO terms contain either the list of the terms or just the
# number of different terms depending on the options.
#
# NOTE: This script does not remove duplicate in the input file. The output file
# should contain rows with UniProt identifiers corresponding to the input file.
#
# Usage:
# python get_GO_terms.py -h



# import
# ------



import argparse



# argument parsing
#-----------------



parser = argparse.ArgumentParser()

parser.add_argument("UniProtIdFile", help = "UniProtIdFile file containing a list of UniProt identifiers",
                    type = str)

parser.add_argument("GOAssociationFile", help = "GOAssociationFile file containing the gene ontology association information",
                    type = str)

parser.add_argument("outputFile", help = "outputFile output file",
                    type = str)

parser.add_argument("-c", "--count", help = "count the number of GO terms in each category instead of returning the list of terms",
                    action = "store_true")

args = parser.parse_args()



# parameters and configuration information
# ----------------------------------------



goa_file = args.GOAssociationFile

UniProt_file = args.UniProtIdFile

output_file = args.outputFile

goa_headers = ["DB", 
               "DB_Object_ID", 
               "DB_Object_Symbol", 
               "Qualifier", 
               "GO_ID", 
               "DB:Reference", 
               "Evidence Code", 
               "With (or) From", 
               "Aspect", 
               "DB_Object_Name", 
               "DB_Object_Synonym", 
               "DB_Object_Type", 
               "Taxon", 
               "Date", 
               "Assigned_By", 
               "Annotation_Extension", 
               "Gene_Product_Form_ID"]

count_GO = args.count



# Functions
#----------



def parseTableToDict(filePath, keys = True, separator = "\t", skipEmpty = True, 
                     returnEmptyOnFileException = False, commentChar = False) :
    
    
    
    """
    
    Takes
    filePath : a path to a file
    keys (True) : a list of keys to be used for the dictionaries
     If set to True, the first line of the file is used as headers
    separator ("\t") : the string used as column separator
    skipEmpty (True) : if True, do not return empty lines in the file 
    commentChar : if set to a character, lines beginning with this character are
      skipped. If a string is given, the first character is used.
    
    Returns
    an iterator over the lines of the file, producing a dictionary 
      per line (non-empty lines only if skipEmpty = True)
    
    """
    
    
    
    fileHandler = open(filePath, "r")
    
    if (keys == True) :
        
        keys = fileHandler.readline().rstrip(" \r\n").split(separator)
        
    for line in fileHandler :
        
        if (not commentChar) :
        
            if (skipEmpty) :
            
                if (line.rstrip() != "") :
            
                    yield(dict(zip(keys, line.rstrip("\r\n").split(separator))))

            else :
                
                yield(dict(zip(keys, line.rstrip("\r\n").split(separator))))
                
        else :
            
            if (not line.startswith(commentChar[0])) :
                
                if (skipEmpty) :
            
                    if (line.rstrip() != "") :
                
                        yield(dict(zip(keys, line.rstrip("\r\n").split(separator))))

                else :
                    
                    yield(dict(zip(keys, line.rstrip("\r\n").split(separator))))

    fileHandler.close()



def loadTableToDictList(filePath, keys = True, keyForKeying = "", separator = "\t", 
                        skipEmpty = True, commentChar = False) :
    
    
    
    """
    
    Load a table from a file and return it as a dictionary using a given field
    from the file as key and a list of dictionariess built from each line as 
    value. This function differs from loadTableToDict since a key with multiple
    entries will have them all in the output, whereas loadTableToDict only return
    the last entry for a given key since multiple entries are overwritting
    each other.
        
    Takes
    filePath : a path to a file
    keys (True) : a list of keys to be used for the dictionaries
     If set to True, the first line of the file is used as headers
    keyForKeying : key of the field to be used as keys for the output dictionary
    separator ("\t") : the string used as column separator
    skipEmpty (True) : if True, do not return empty lines in the file 
    commentChar : if set to a character, lines beginning with this character are
      skipped. If a string is given, the first character is used.

    Returns
    a dictionary containing the table
    
    """
    
    
    
    fileIterator = parseTableToDict(filePath = filePath, keys = keys, 
                                    separator = separator, skipEmpty = skipEmpty,
                                    commentChar = commentChar)
    
    output = {}
    
    for line in fileIterator :
        
        entry = output.get(line[keyForKeying], [])
        
        entry.append(line)
        
        output[line[keyForKeying]] = entry
        
    return output






# script
#-------



# load the GO table

GO_dict = loadTableToDictList(goa_file, keys = goa_headers,
            keyForKeying = "DB_Object_ID", commentChar = "!")



# load the UniProt identifier list

fi = open(UniProt_file, "r")

UniProt_id = []

for line in fi :
    
    if (line.strip() != "") :
        
        UniProt_id.append(line.strip())



# write the full GO term table for UniProt id

f = open(output_file, "w")

# prepare the headers

f_headers = ["UniProtId", "GO_BP", "GO_CC", "GO_MF"]

f.write("\t".join(f_headers) + "\n")

# go through the UniProt identifier list

for UniProt_each_id in UniProt_id :

    # get the GO information for that UniProt identifier

    GO_list = GO_dict.get(UniProt_each_id, [])
    
    # prepare the empty lists for each aspect (P, C, F)
    
    aspect_dict = {}
    
    aspect_dict["P"] = []
    
    aspect_dict["C"] = []
    
    aspect_dict["F"] = []
    
    # assign each GO term to the corresponding list
    
    for GO_entry in GO_list :
        
        aspect_dict[GO_entry["Aspect"]].append(GO_entry["GO_ID"])
    
    # remove duplicates in GO_ID
    
    for aspect in ["P", "C", "F"] :
        
        aspect_dict[aspect] = list(set(aspect_dict[aspect]))
    
    # write it to the file
    
    if (not count_GO) :
    
        # write the detailed list of terms
    
        output = [UniProt_each_id, 
                  ",".join(aspect_dict["P"]), 
                  ",".join(aspect_dict["C"]),
                  ",".join(aspect_dict["F"])]
                  
        f.write("\t".join(output) + "\n")

    else :
        
        # count the number of terms
        
        output = [UniProt_each_id, 
                  str(len(aspect_dict["P"])), 
                  str(len(aspect_dict["C"])),
                  str(len(aspect_dict["F"]))]
                  
        f.write("\t".join(output) + "\n")
    
f.close()
