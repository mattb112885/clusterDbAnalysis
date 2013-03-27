#!/usr/bin/python

# Get gene neighborhoods within +/- K (default = 3) from each gene in stdin
#

import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *

usage = "%prog [options] < gene_id_list > gene_neighborhoods"
description="""Given a list of gene IDs, get the neighborhoods within the specified 
number of genes on the same contig on either strand from the specified gene.
Note that start location is always the first base of a start codon (so start > stop for - strand genes).
The neighborhoods table has the following fields:
querygene:::neighbor_gene:::gene_distance:::contig:::start_location:::strand:::annotation"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--neighborhood", help="Size of desired neighborhood in number of genes from center (maximum of 5, D=3)", action="store", type="int", dest="nsize", default=3)
parser.add_option("-g", "--genecol", help="Column number starting from 1 for gene ID (D=1)", action="store", type="int", dest="gc", default=1)
(options, args) = parser.parse_args()

gc = options.gc - 1

nsize = options.nsize
if nsize > 10:
    sys.stderr.write("ERROR: Maximum neighborhood size is 10 (this is the size used when pre-calculating for input into the database\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for gene in fileinput.input("-"):
    spl = gene.strip('\r\n').split("\t")
    geneid = spl[gc]
    if not geneid.startswith("fig|"):
        geneid = "fig|%s" %(geneid)
    cur.execute("""SELECT * from neighborhoods
                   WHERE neighborhoods.centergene=? AND ABS(neighborhoods.distance) <= ?;""", (geneid,nsize))
    for l in cur:
        print "\t".join([ str(s) for s in list(l) ])

con.close()
