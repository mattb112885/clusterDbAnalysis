#!/usr/bin/env python

import optparse, sys, sqlite3
from FileLocator import *

ok_databases = ["all", "cd", "cog", "pfam", "tigr", "prk", "smart"]
header = [ "cdd_id", "external_clusterid", "clustername", "description", "profile_length" ]

usage = """%prog [Options] \"Description 1\" \"Description 2\" ...

Output table:  """ + " ".join(header)
description = """Identify external clusters (optionally in a specific database)
that match a specified description. If more than one description is provided, return
all clusters matching at least one description."""
parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-d", "--database", help="External clustering database to search. Default: all those in NCBI CDD. Valid options are: %s" %(" ".join(ok_databases)),
                  action="store", type="str", dest="database", default="all")
(options, args) = parser.parse_args()

if len(args) == 0:
    sys.stderr.write("ERROR: At least one description to match is required (use -h for usage details) \n")
    exit(2)

db = options.database.lower()

if db not in ok_databases:
    sys.stderr.write("ERROR: Invalid database passed with -d flag. Valid databases are: %s\n" %(" ".join(ok_databases)))
    exit(2)

# Set up SQL query inputs...
teststr = list('%' + s + '%' for s in args)

# Set up SQL query
query = """SELECT * FROM external_clusters WHERE """
if db != "all":
    query += " external_clusterid LIKE ? AND "

toAddToQuery = []
for ii in range(len(teststr)):
    toAddToQuery.append(" clustername LIKE ? OR description LIKE ? OR external_clusterid LIKE ? ")
toAddToQuery = "(" + "OR".join(toAddToQuery) + ")"
query += toAddToQuery

toSubstitute = []
if db != "all":
    toSubstitute.append('%' + db + '%')
for st in teststr:
    # This isn't a typo - we need to add the same one thrice in a row
    toSubstitute.append(st)
    toSubstitute.append(st)
    toSubstitute.append(st)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()
cur.execute(query, tuple(toSubstitute))

for rec in cur:
    print "\t".join(str(s) for s in rec)
