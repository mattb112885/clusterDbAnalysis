#!/usr/bin/env python

from __future__ import print_function
import fileinput
import optparse
import sqlite3
import sys

from ClusterGraph import *
from FileLocator import *

symmetric_methods, unsymmetric_methods = getValidBlastScoreMethods()
all_methods = symmetric_methods + unsymmetric_methods

usage = """ %prog -m method -u cutoff [options] < BLAST_results """
description = """Generate a GML file for all edges in a BLAST results table.
Currently implemented methods: %s """ %(" ".join(all_methods) )

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--method", help="Scoring method to use (REQUIRED)", action="store", type="str", dest="method", default=None)
parser.add_option("-u", "--cutoff", help="Scoring cutoff to use (REQUIRED)", action="store", type="float", dest="cutoff", default=None)
parser.add_option("-o", "--outfile", help="Output GML file (REQUIRED)", action="store", type="str", dest="outfile", default=None)
(options, args) = parser.parse_args()

if options.method is None or options.cutoff is None or options.outfile is None:
    sys.stderr.write("ERROR: method (-m), cutoff (-u) and output file (-o) are required arguments\n")
    exit(2)

blastres = [ line.strip("\r\n").split("\t") for line in fileinput.input("-") ]

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

G = makeNetworkObjectFromBlastResults( blastres, options.method, options.cutoff, cur )

exportGraphToGML(G, options.outfile)
