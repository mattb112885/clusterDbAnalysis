#!/usr/bin/env python

import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *
from ClusterFuncs import *

header = [ "contig", "nucleotide_sequence" ]
usage = """%prog [options] < contig_ids > contig_sequences

Output : """ + " ".join(header) + "\n" + """Alternative output: FASTA file"""
description = """Get the complete DNA sequence for contigs with specified IDs."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--fasta", help="Return sequences as FASTA file (the default is a two-column table)", 
                  action="store_true", dest="fasta", default=False)
parser.add_option("-c", "--contigcol", help="Column number for contig ID, starting from 1 (D;1)",
                  action="store", dest="cc", default=1)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

cc = options.cc - 1

contiglist = set()
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    contiglist.add(spl[cc])

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

contigdict = getContigSequence(cur, contiglist)
if options.header:
    print "\t".join(header)
for contig in contigdict:
    if options.fasta:
        print ">%s\n%s\n" %(contig, contigdict[contig])
    else:
        print "%s\t%s" %(contig, contigdict[contig])

con.close()
