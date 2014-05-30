#!/usr/bin/env python

# Provide a list of organisms to match [can match any portion of the organism
# so if you give it just "mazei" it will return to you a list of Methanosarcina mazei]
#
# Returns a list of blast results specific to those organisms to stdout
# (which can subsequently be used to do clustering...)
# 

import sqlite3, optparse, fileinput
from FileLocator import *
from ClusterFuncs import *

header = [ 'Query_gene', 'Target_gene', 'Percent_identity', 'HSP_length', 'Percent_mismatch',
           'Gap_opens', 'Query_start', 'Query_end', 'Target_start', 'Target_end', 'E_value',
           'Bit_score', 'Query_selfbit', 'Target_selfbit' ]

usage = """%prog [options] < gene_ids > blast_results

Output: """ + " ".join(header)

description = "Given list of genes to match, returns a list of BLAST results between genes in the list only"

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
parser.add_option("-n", "--blastn", help="Get BLASTN results (D: BLASTP results)", action="store_true", dest="blastn", default=False)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()
gc = options.genecolumn - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

geneids = []
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    geneids.append(spl[gc])

blastres = getBlastResultsBetweenSpecificGenes(geneids, cur, blastn=options.blastn)
if options.header:
    print "\t".join(header)
for res in blastres:
    print "\t".join(res)

con.close()
