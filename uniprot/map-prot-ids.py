#! /usr/bin/python

Version = "0.0.1"

### * Setup

### ** Import

import argparse
import math
import logging

### ** Logger

log = logging.getLogger("log")
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)

### ** Parameters

maxRecordsPerBatch = 5000 # set arbitrarily, nothing found about it on uniprot.org
allowed_db = ["ACC", "ID", "UPARC", "NF50", "NF90", "NF100",
              "EMBL_ID", "EMBL", "PIR", "UNIGENE_ID", "P_ENTREZGENEID",
              "P_GI", "P_REFSEQ_AC", "REFSEQ_NT_ID"]

### * Argument parsing

### ** Description

description = ("GI-to-UniProtKB.py version " + Version + "\n" +
  "Convert GI identifiers to UniProtKB identifiers, using code from Jan Rudolph")
epilog = ("Accepted abbrevations are (see http://www.uniprot.org/help/programmatic_access):" + "\n" +
          " -Uniprot-\n" +
          "  ACC...............UniProtKB AC\n" +
          "  ID................UniProtKB ID\n" +
          "  UPARC.............UniParc\n" +
          "  NF50..............UniRef50\n" +
          "  NF90..............UniRef90\n" +
          "  NF100.............UniRef100\n" +
          " -Other sequence databases-\n" +
          "  EMBL_ID...........EMBL/GenBank/DDBJ\n" +
          "  EMBL..............EMBL/GenBank/DDBJ CDS\n" +
          "  PIR...............PIR\n" +
          "  UNIGENE_ID........UniGene\n" +
          "  P_ENTREZGENEID....Entrez Gene (GeneID)\n" +
          "  P_GI..............GI number\n" +
          "  P_REFSEQ_AC.......RefSeq Protein\n" +
          "  REFSEQ_NT_ID......RefSeq Nucleotide")

### ** Parser

parser = argparse.ArgumentParser(description = description, epilog = epilog,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", "--input", metavar = "INPUT_FILE",
                    help = "input file (one column, no header)",
                    type = str)
parser.add_argument("-o", "--output", metavar = "OUTPUT_FILE",
                    help = "output file (two columns, no headers)",
                    type = str)
parser.add_argument("-b", "--batchSize", metavar = "N",
                    help = "number of records per batch submitted to uniprot.org",
                    type = int, default = 5000)
parser.add_argument("-f", "--db_from", metavar = "DB_ABBR",
                    help = "abbreviation of the database to map from",
                    type = str)
parser.add_argument("-t", "--db_to", metavar = "DB_ABBR",
                    help = "abbreviation of the database to map to",
                    type = str)
parser.add_argument("--email", metavar = "USER_EMAIL",
                    help = "user's email, provided to NCBI services",
                    type = str)


### ** Parse the arguments

args = parser.parse_args()
if (args.batchSize > maxRecordsPerBatch) :
    raise Exception("batchSize is too high, max value is " + str(maxRecordsPerBatch))
if (not args.email) :
    raise Exception(("You should provide an email address to be sent to UniProt " +
                     "with the requests"))
recordsPerBatch = args.batchSize
contact_email = args.email
db_from = args.db_from
db_to = args.db_to
assert db_from in allowed_db
assert db_to in allowed_db

### * Third-party code

"""
code from Jan Rudolph
from https://pypi.python.org/pypi/uniprot_tools/0.4.1

modified by Matthieu Bruneaux

uniprot python interface
to access the uniprot database

available services:
    map
    retrieve
"""

### ** Setup

import requests
import sys

url = 'http://www.uniprot.org/'

### ** Functions

### *** _retrieve(query, format='txt')

def _retrieve(query, format='txt'):
    """_retrieve is not meant for use with the python interface, use `retrieve`
    instead"""
    tool = 'batch/'

    query = list(set(query.split()))
    queries = [query[i:i+100] for i in xrange(0, len(query), 100)]

    data = {'format':format}

    responses = [requests.post(url + tool,
                               data = data,
                               files = {'file':' '.join(query)},
                               header = {"From" : contact_email}) for query in queries]
    page = [response.text for response in responses]
    return page

### *** retrieve(ids, format='txt')

def retrieve(ids, format='txt'):
    """ request entries by uniprot acc using batch retrieval

    Args:
        query: list of ids to retrieve
        format: txt by default

    Help:
        possible formats:
        txt, xml, rdf, fasta, gff"""
    if type(ids) is not list:
        ids = [ids]
    return _retrieve(' '.join(ids), format)

### *** _map(query, f, t, format='tab')

def _map(query, f, t, format='tab'):
    """ _map is not meant for use with the python interface, use `map` instead
    """
    tool = 'mapping/'

    data = {
            'from':f,
            'to':t,
            'format':format,
            'query':query
            }
    response = requests.post(url + tool,
                             data = data,
                             headers = {"From" : contact_email})
    page = response.text
    return page

### *** map(ids, f, t, format='tab')

def map(ids, f, t, format='tab'):
    """ map a list of ids from one format onto another using uniprots mapping api
    
    Args:
        query: id or list of ids to be mapped
        f: from ACC | P_ENTREZGENEID | ...
        t: to ...
        format: tab by default

    Help:
        for a list of all possible mappings visit
        'http://www.uniprot.org/faq/28'
    """
    if type(ids) is not list:
        ids = [ids]
    page = _map(' '.join(ids), f, t, format)
    result = dict()
    for row in page.splitlines()[1:]:
        key, value = row.split('\t')
        if key in result:
            result[key].add(value)
        else:
            result[key] = set([value])
    return result

### * Functions

### ** loadGI(inputFile)

def loadGI(inputFile) :
    with open(inputFile, "r") as f :
        GIs = []
        for l in f :
            l = l.strip()
            if (l != "") :
                GIs.append(l)
        return(GIs)
    raise Exception("Error while reading input file: " + fileName)

### * Run

### ** Load the GI identifiers

log.info("Loading the GI identifiers")
GI_ids = list(set(loadGI(inputFile = args.input)))
log.info("There are " + str(len(GI_ids)) + " unique GI.")

### ** Create the batches of GI to submit

log.info("Creating the submission batches")
nGIs = len(GI_ids)
nBatches = int(math.ceil(nGIs * 1. / recordsPerBatch))
batches = [GI_ids[(i * recordsPerBatch) : min(((i + 1) * recordsPerBatch), nGIs)]
           for i in range(nBatches)]
# check
rGIs = []
for x in batches :
    rGIs += x
assert rGIs == GI_ids

### ** Submit the mapping to UniProt services

log.info("Mapping the GI identifiers to UniProtKB identifiers using uniprot.org")
log.info("Mapping from " + db_from + " to " + db_to)
mapping = dict()
for i in range(nBatches) :
    log.info("Requesting batch number " + str(i + 1) + "/" + str(nBatches))
    r = map(batches[i], f = db_from, t = db_to, format = "tab")
    mapping.update(r)

### ** Write the results

log.info("Writing the results")
with open(args.output, "w") as f :
    for (k, v) in mapping.items() :
        l = "\t".join([k] + list(v)) + "\n"
        f.write(l)

### * End

log.info("Complete")
