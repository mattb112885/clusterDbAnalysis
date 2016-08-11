#!/usr/bin/env python

import fileinput
import optparse
import re
import sqlite3
import sys

from FileLocator import *
from sanitizeString import *
from ClusterFuncs import *

usage = """%prog [options] focus_gene < gene_ids > cluster_comparison_table"""
description = '''
Given a list of gene IDs and a focus gene, identifies all of the clusters (across
all cluster runs in the database) that contain the gene. Then identifies if the
genes in 'gene_ids' are identified in the same cluster as the focus gene or not. 
The output is a table with 1 if the gene is present and 0 if it is not (some options 
allow changing this to yes/no). An example usage is to provide all genes in a broad cluster (low
cutoff) and use this to see how cutoff affects the clustering results. The evaluation is done for every
cluster run in which both genes are found.

Provided gene IDs can be sanitized or unsanitized.'''

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-n", "--newick", help="Specify this if input is a Newick file and not a list of gene IDs", action="store_true", dest="newick", default=False)
parser.add_option("-g", "--genecol", help="Column number for gene ID in input starting from 1 (not valid with -n, D=1)", action="store", dest="gc", type="int",  default=1)
parser.add_option("-y", "--yesno", help="Print yes and no instead of 1 and 0", action="store_true", dest="yesno", default=False)
parser.add_option("-f", "--filter_runs", help="Only include cluster runs matching the specified regex", action="store", dest="filter_runs", default=None)
parser.add_option("-a", "--all", help="Instead of specifying a list of gene IDs, this will print a table for ALL genes found in the same cluster as the focus gene in at least one run.",
                  action="store_true", dest="all", default=False)
(options, args) = parser.parse_args()

if len(args) < 1:
    sys.stderr.write("ERROR: Focus gene is a required input\n")
    exit(2)

if options.filter_runs is not None:
    options.filter_runs = re.compile(options.filter_runs)

focusgene = args[0]

gc = options.gc - 1

# Get a list of cluster/runID pairs containing the focus gene
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
res = getClustersContainingGenes( [ focusgene ], cur)
# Just in case there is more than one cluster per run
clusterrun_to_genes = {}
for result in res:
    runid = result[0]
    clusterid = result[1]
    if options.filter_runs is not None and options.filter_runs.search(runid) is None:
        continue        
    geneids = getGenesInCluster(runid, clusterid, cur)
    clusterrun_to_genes[(runid, clusterid)] = set(geneids)

# Read in the remainder of the gene IDs
geneids = set()
if options.all:
    geneids = []
    for clusterrun in clusterrun_to_genes:
        geneids += list(clusterrun_to_genes[clusterrun])
    geneids = set(geneids)
elif options.newick:
    from ete2 import Tree
    from TreeFuncs import *
    newick_str = "".join( [ line.strip("\r\n") for line in fileinput.input("-") ] )
    t = Tree(newick_str)
    geneids = set(getLeafNames(t))
else:
    for line in fileinput.input("-"):
        spl = line.strip("\r\n").split("\t")
        geneids.add(spl[gc])

# Set up the table. We set it up to hopefully plug right in to the
# Newick file if that is what was provided
header = ['#names']
for clusterrun in sorted(clusterrun_to_genes.keys()):
    cr_id = "%s_%s" %(clusterrun[0], clusterrun[1])
    header.append(cr_id)
print("\t".join(header)) 

# Now we look for each gene in each cluster.
for gene in geneids:
    # In order to work with db_displayTree this should match with what is
    # in the Newick tree (whether it is sanitized or not)
    row = [ gene ]
    gene = unsanitizeGeneId(gene)

    for clusterrun in sorted(clusterrun_to_genes.keys()):
        if gene == focusgene and options.yesno:
            printstr = "focus"
        elif gene in clusterrun_to_genes[clusterrun]:
            if options.yesno:
                printstr = "yes"
            else:
                printstr = "1"
        else:
            if options.yesno:
                printstr = "no"
            else:
                printstr = "0"
        row.append(printstr)
    print("\t".join(row))

con.close()
