#!/usr/bin/python

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
# echo "fig|188937.5.peg.1" | ./src/db_getGeneInformation
# will give you all the info on gene fig|188937.5.peg.1
#
# If you do not include the "fig|" I will add it for you. However, no other modifications are allowed.
#


import fileinput, sqlite3, optparse

usage="%prog [options] < gene_ids > gene_info"
description="Given a list of gene IDs, get their gene info, including annotations, contig, organism, strand, and sequences"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
(options, args) = parser.parse_args()
gc = options.genecolumn - 1 # Convert to Pythonic indexes                                                                                 

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip().split("\t")

    # Allow shortcuts if we start by looking at the seedviewer website for the fig ID
    if not spl[gc].startswith("fig|"):
        spl[gc] = "fig|" + spl[gc]

    cur.execute("SELECT * FROM processed WHERE processed.geneid = ?;", (spl[gc],))

    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri
    
con.close()
