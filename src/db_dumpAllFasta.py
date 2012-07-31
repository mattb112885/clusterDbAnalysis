#!/usr/bin/python

# Dump a fasta file with UPDATED annotations (including ALIASES)
# from the database.
#
# The fasta file will contain all of the genes in the database...

import optparse, sqlite3, sys
from locateDatabase import *

usage="$prog > Fasta_file"
description="Generates a fasta file with all the annotations in the database including aliases added to the raw annotations..."
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--nucleotide", help="Make nucleotide fasta file (D: Protein)", action="store_true", dest="nt", default=False)
(options, args) = parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

if options.nt:
    query = "SELECT geneid,organism,annotation,nucseq FROM processed;"
else:
    query = "SELECT geneid,organism,annotation,aaseq FROM processed;"

cur.execute(query)
for res in cur:
    ls = [ str(s) for s in res ]
    print ">%s_%s_%s\n%s" %(ls[0], ls[1], ls[2], ls[3])

con.close()
