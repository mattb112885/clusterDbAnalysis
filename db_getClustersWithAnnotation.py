#!/usr/bin/python

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

import sqlite3, fileinput, optparse

usage = "%prog \"Annotation 1\" \"Annotation 2\" ..."
description = "Given list of run IDs, returns a list of genes and clusters containing given word(s) in the annotation"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

# Change the annotations so that they are all LIKE %name%
teststr = tuple('%' + s + '%' for s in args)

con = sqlite3.connect("db/methanosarcina")
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

for line in fileinput.input("-"):
    # Note - "+" for tuples is the concatination operator
    cur.execute(query, (line.strip(),) + teststr)
    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri

con.close()
fileinput.close()
