#!/usr/bin/python

# Get gene neighborhoods within +/- K (default = 3) from each gene in stdin
#

import fileinput
import sqlite3
import optparse
import sys
from locateDatabase import *

usage = "%prog [options] < gene_id_list > gene_neighborhoods"
description="Given a list of gene IDs, get the neighborhoods within the specified number of genes on the same contig on either strand from the specified gene"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--neighborhood", help="Size of desired neighborhood in number of genes from center (maximum of 5, D=3)", action="store", type="int", dest="nsize", default=3)
parser.add_option("-g", "--genecol", help="Column number starting from 1 for gene ID (D=1)", action="store", type="int", dest="gc", default=1)
(options, args) = parser.parse_args()

gc = options.gc - 1

nsize = options.nsize
if nsize > 5:
    sys.stderr.write("ERROR: Maximum neighborhood size is 5 (this is the size used when pre-calculating for input into the database\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for gene in fileinput.input("-"):
    spl = gene.strip('\r\n').split("\t")
    cur.execute("""SELECT neighborhoods.*, processed.annotation FROM neighborhoods
                   INNER JOIN processed ON processed.geneid = neighborhoods.neighborgene
                   WHERE neighborhoods.centergene=? AND ABS(neighborhoods.distance) <= ?;""", (spl[gc],nsize))
    for l in cur:
        print "\t".join([ str(s) for s in list(l) ])

con.close()
