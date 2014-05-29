#!/usr/bin/env python

import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *
from ClusterFuncs import *

usage = """
%prog runID -n 'organism_name' < gene_list > equivalent_genes
%prog runID -i 'organism_id'   < gene_list > equivalent_genes
"""

description = """Given a list of genes and an organism name or ID, identifies all genes
in the same cluster as the target genes in the specified organism.

Returns a two-column table containing the list of genes in the first column
and the equivalent gene IDs in the other organism separated by semicolons in the second. 
Note that the order of the output is not necessarily the same as the order of the input.
"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--name", help="Name of target organism.", action="store", type="str", dest="orgname", default=None)
parser.add_option("-i", "--id", help="Organism ID for target organism.", action="store", type="str", dest="orgid", default=None)
parser.add_option("-g", "--genecol", help="Column number for gene Id starting from 1 (D:1)", action="store", type="int", dest="gc", default=1)
(options, args) = parser.parse_args()

if len(args) < 1:
    sys.stderr.write("ERROR: run ID is a required argument.\n")
    exit(2)

runid = args[0]

gc = options.gc - 1

genelist = set()
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    genelist.add(spl[gc])

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

equivalent_dict = getEquivalentGenesInOrganism( genelist, runid, cur, orgid=options.orgid, orgname=options.orgname )

for gene in equivalent_dict:
    print "\t".join( [ gene, ";".join(equivalent_dict[gene]) ] )

con.close()
