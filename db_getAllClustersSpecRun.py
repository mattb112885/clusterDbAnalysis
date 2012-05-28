#!/usr/bin/python

# This is a pipe command
#
# Pipe in a run ID (or list of run IDs) and it returns to you a 2-column table containing
# the run ID(s) in the first column and all of the cluster IDs contained within it in the second
# column (and gene IDs in the third column)
#
# Results are printed to stdout.

import fileinput, sqlite3, optparse, os

usage="%prog"
description="Given a set of run IDs (from stdin), returns all cluster IDs associated with that run ID"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

for line in fileinput.input("-"):
    runid = line.strip('\n')
    query = "SELECT * FROM clusters WHERE clusters.runid = ?;"
    cur.execute(query, (runid, ) )
    
    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri

con.close()
