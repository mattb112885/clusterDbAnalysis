#!/usr/bin/env python

# This is a pipe command
#
# Pipe in a stdin with genes located either in the 1st column (default)
# or in the column deonted by -g
#
# It will output all available information (from the PROCESSED table) for those genes 
# by default.
#
# You need to use the entire ID, i.e.
#
# echo "fig|188937.5.peg.1" | db_getGeneInformation.py
# will give you all the info on gene fig|188937.1.peg.1
#
# If you do not include the "fig|" I will add it for you. However, no other modifications are allowed.
#

from __future__ import print_function
import fileinput, sqlite3, optparse
from FileLocator import *

headers = [ "organism_name", "contig_id", "start", "stop", "strand", "strandnum", "annotation", "DNA_seq", "AA_seq" ]

usage="""%prog [options] < gene_ids > gene_info

Output table: """ + " ".join(headers)
description="""Given a list of gene IDs, get their gene info, 
including annotations, contig, organism, strand, and sequences. 
Start is the first nucleotide of the start codon (e.g. A in ATG)
and stop is the last nucleotide, relative to nucleotide 1 at the start
of the contig.
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
parser.add_option("-a", "--add", help="Add gene information to the end of the existing file (D: only return the gene information)", 
                  action="store_true", dest="keep", default=False)
(options, args) = parser.parse_args()

gc = options.genecolumn - 1 # Convert to Pythonic indexes               

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")

    # Allow shortcuts if we start by looking at the seedviewer website for the fig ID
    if not spl[gc].startswith("fig|"):
        spl[gc] = "fig|" + spl[gc]

    cur.execute("SELECT processed.* FROM processed WHERE processed.geneid = ?;", (spl[gc],))

    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        if options.keep:
            stri = "%s\t%s" %(line.strip('\r\n'), stri)
        print(stri)
    
con.close()
