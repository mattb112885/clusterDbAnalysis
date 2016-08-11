#!/usr/bin/env python

# This is a pipe command.
#
# Unlike db_getBlastResultsBetweenSpecificGenes (which only gets blast results containing both
# query and target inside them),
# this command will get blast results if either the query OR the target
# is contained within them.
#
# The code is identical except the sql command has an OR instead of an AND...

import fileinput, optparse, sqlite3
from ClusterFuncs import *
from FileLocator import *

usage = "%prog [options] < gene_ids > blast_results"
description = "Given list of genes to match, returns a list of BLAST results containing any gene ID in your list as either a query or a target (for blast results only BETWEEN the query genes, see db_getBlastResultsBetweenSpecificGenes.py)"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
parser.add_option("-c", "--cutoff", help="E-value cutoff (D: Show all results in database)", action="store", type="float", dest="cutoff", default=10)
parser.add_option("-n", "--blastn", help="Base the results on BLASTN instead of BLASTP (D: BLASTP)", action="store_true", dest="blastn", default=False)
parser.add_option("-o", "--onlyquery", help="Only return results for which one of the specified genes is the query (by default returns results whether it is a query or a target", action="store_true",
                  dest="only_query", default=False)
(options, args) = parser.parse_args()

gc = options.genecolumn - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

genelist = set()
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    gn = spl[gc]
    if not gn.startswith("fig|"):
        gn = "fig|%s" %(gn)
    genelist.add(gn)

blastres = getBlastResultsContainingGenes(genelist, cur, cutoff=options.cutoff, blastn=options.blastn, only_query=options.only_query)

for blast in blastres:
    print("\t".join(blast))

con.close()
