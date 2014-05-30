#!/usr/bin/env python

# This is a pipe command
#
# Pipe in run ID(s) (from db_getClusterId)
#
# Generates a table
# runid | clusterid | geneid | annotation
#
# For example:
# ./db_getClustersWithAnnotation.py "Hypothetical"
# 
# would list all genes with "Hypothetical" in their names and what cluster
# they belonged to (note - it does NOT return all of the genes belonging to each of these clusters)
#
# The matching is not case sensitive and you can provide as many annotations to match as you want.
#
# Multiple different annotations are combined with OR statements.

import sqlite3, fileinput, optparse, sys
from FileLocator import *

header = [ "Runid", "Cluster_id", "Matching_gene", "Matching_annotation" ]
usage = """%prog [options] \"Annotation 1\" \"Annotation 2\" ... < run_ids > clusters_with_genes_containing_annotation_words

Output: """ + " ".join(header)
description = """Given list of run IDs, returns a list of genes and clusters containing at least one of the given word(s) in the annotation.
Note that only the genes matching the annotation are returned, not all of the genes in the clusters. Use db_getGenesInClusters.py to get
all of the genes in a list of clusters."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--runcol", help="Column number for run ID, starting from 1 (D=1)", action="store", type="int", dest="rc", default=1)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

if len(args) < 1:
    sys.stderr.write("ERROR: Must provide at least one annotation to match!\n")
    exit(2)

rc = options.rc - 1

# Change the annotations so that they are all LIKE %name%
teststr = tuple('%' + s + '%' for s in args)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

query = """SELECT clusters.*, processed.annotation
           FROM clusters
           INNER JOIN processed ON clusters.geneid = processed.geneid
           WHERE clusters.runid = ? AND 
           ("""

for i in range(len(teststr)):
    query = query + "processed.annotation LIKE ? "
    if not i == len(teststr) - 1:
        query = query + "OR "
query = query + ");"

if options.header:
    print "\t".join(header)
for line in fileinput.input("-"):
    # Note - "+" for tuples is the concatination operator
    spl = line.strip('\r\n').split("\t")
    cur.execute(query, (spl[rc],) + teststr)
    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri

con.close()
fileinput.close()
