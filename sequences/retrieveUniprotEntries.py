#! /usr/bin/python

# Description
# -----------
#
# This script enables to perform some id mapping using the UniProt services
# (see http://www.uniprot.org/uploadlists/ and
# http://www.uniprot.org/help/programmatic_access).

Version = "0.0.1"

### * Setup

### ** Imports

import argparse
import urllib
import urllib2

### ** Parameters

url_uniprot = "http://www.uniprot.org/mapping/"
user_email = "mdjbru@utu.fi"

### * Argument parsing

### ** Description

description = ("mapUniProtId.py version " + Version + "\n" +
               "\n" +
               "Valid abbreviations for the databases can be found at \n" +
               "http://www.uniprot.org/help/programmatic_access#id_mapping_examples \n" +
               "The script was modified from the example Python script available on this " +
               "page.")

### ** Parser

parser = argparse.ArgumentParser(description = description)
parser.add_argument("-f", "--db_from", metavar = "DATABASE_FROM",
                    help = "abbreviation of the database name to convert from",
                    type = str)
parser.add_argument("-t", "--db_to", metavar = "DATABASE_TO",
                    help = "abbreviation of the database name to convert to",
                    type = str)
parser.add_argument("ids", metavar = "IDENTIFIER", nargs = "+",
                    help = "one or several identifiers to convert",
                    type = str)

### ** Parse the arguments

args = parser.parse_args()

### * Run

url_parameters = {"from" : args.db_from,
                  "to" : args.db_to,
                  "format" : "tab",
                  "query" : " ".join(args.ids)}
url_data = urllib.urlencode(url_parameters)
url_request = urllib2.Request(url_uniprot, url_data)
url_contact = user_email
url_request.add_header("User-Agent", "Python %s" % url_contact)
response = urllib2.urlopen(url_request)
page = response.read(200000)

f = open("toto", "w")
f.write(page)
f.close()

### * Test

# python mapUniProtId.py -f P_GI -t ID 443713586
