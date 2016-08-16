#!/usr/bin/env python

from __future__ import print_function
import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *
from ClusterFuncs import *

usage = "%prog [options] < contig_ids > contig_sequences"
description = """Get the complete DNA sequence for contigs with specified IDs."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--fasta", help="Return sequences as FASTA file (the default is a two-column table)", 
                  action="store_true", dest="fasta", default=False)
parser.add_option("-c", "--contigcol", help="Column number for contig ID, starting from 1 (D;1)",
                  action="store", dest="cc", default=1)
(options, args) = parser.parse_args()

cc = options.cc - 1

contiglist = set()
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    contiglist.add(spl[cc])

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

contigdict = getContigSequence(cur, contiglist)
for contig in contigdict:
    if options.fasta:
        print(">%s\n%s\n" %(contig, contigdict[contig]))
    else:
        print("%s\t%s" %(contig, contigdict[contig]))

con.close()
