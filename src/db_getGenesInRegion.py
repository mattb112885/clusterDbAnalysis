#!/usr/bin/python

import fileinput
import optparse
import sqlite3

from ClusterFuncs import *
from FileLocator import *

usage = "%prog [options] < contig_and_gene_locations > geneinfo"
description = """ This function takes a table with contig and regions (start and stop)
as inputs and identifies all of the genes between the start and stop columns (by default, the entire 
gene must be between the start and stop). It returns the geneinfo for the genes within
the specified regions. The defaults are set up to take a TBLASTn results file from the ITEP
TBLASTn wrapper as input"""

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-c", "--contigcol", help="Column number for the contig to search, starting from 1 (D=3)",
                  action="store", type="int", dest="contigcol", default=3)
parser.add_option("-s", "--startcol", help="Column number for the start of the region to search, starting from 1 (D=5)",
                  action="store", type="int", dest="startcol", default=5)
parser.add_option("-e", "--endcol", help="Column number for the end of the region to search, starting from 1 (D=6)",
                  action="store", type="int", dest="endcol", default=6)
parser.add_option("-x", "--expand", 
                  help="Expand the gene regions by this many nucleotides before searching (D: 0 - use regions as provided)",
                  action="store", type="int", dest="expand", default=0)
parser.add_option("-v", "--overhang",
                  help = "Allow this number of amino acids to hang over the edge of the expanded region by this many nucleotides (D: 0 - the entire gene must lay in the region)",
                  action = "store", type="int", dest="overhang", default=0)

(options, args) = parser.parse_args()

# Convert to python indexes
cc = options.contigcol - 1
sc = options.startcol - 1
ec = options.endcol - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    contig = spl[cc]
    start = min(int(spl[sc]), int(spl[ec])) - options.expand
    stop = max(int(spl[sc]), int(spl[ec])) + options.expand

    geneids = getGenesInRegion(contig, start, stop, cur, overhang=options.overhang)
    geneinfo = getGeneInfo(geneids, cur)
    for info in geneinfo:
        print "\t".join(info)

con.close()
