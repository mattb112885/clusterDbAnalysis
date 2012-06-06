#!/usr/bin/python

# Get gene neighborhoods within +/- K (default = 3) from each gene in stdin
#

import fileinput
import sqlite3
import optparse
import sys

usage = "%prog [options] < gene_id_list > gene_neighborhoods"
description="Given a list of gene IDs, get the neighborhoods within the specified number of genes on the same contig on either strand from the specified gene"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--neighborhood", help="Size of desired neighborhood in number of genes from center (maximum of 5)", action="store", type="int", dest="nsize", default=3)
(options, args) = parser.parse_args()

nsize = options.nsize
if nsize > 5:
    sys.stderr.write("ERROR: Maximum neighborhood size is 5 (this is the size used when pre-calculating for input into the database\n")
    exit(2)

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

for gene in fileinput.input("-"):
    cur.execute("""SELECT neighborhoods.*, processed.annotation FROM neighborhoods
                   INNER JOIN processed ON processed.geneid = neighborhoods.neighborgene
                   WHERE neighborhoods.centergene=? AND ABS(neighborhoods.distance) <= ?;""", (gene.strip(),nsize))
    for l in cur:
        print "\t".join([ str(s) for s in list(l) ])

con.close()
