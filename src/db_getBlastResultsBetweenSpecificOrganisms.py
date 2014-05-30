#!/usr/bin/env python

# Provide a list of organisms to match [can match any portion of the organism
# so if you give it just "mazei" it will return to you a list of Methanosarcina mazei]
#
# This function returns only those that match ALL of the criteria! To get all of the organisms put some nonsense like
# !ZZZ
#
# which returns all that do NOT match ZZZ (which is all of them presumably!)
#
# Returns a list of blast results specific to those organisms to stdout
# (which can subsequently be used to do clustering...)
# 

import sqlite3, optparse
from FileLocator import *

header = [ 'Query_gene', 'Target_gene', 'Percent_identity', 'HSP_length', 'Percent_mismatch',
           'Gap_opens', 'Query_start', 'Query_end', 'Target_start', 'Target_end', 'E_value',
           'Bit_score', 'Query_selfbit', 'Target_selfbit' ]


usage = """%prog \"Organism 1\" \"Organism 2\" ... > blast_results

Output: """ + " ".join(header)

description = """Given list of organism names to match, returns a list
of BLAST results between organisms matching any of those keywords."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-s", "--strict", help="Require exact name matches for organisms (D: Partial name matches are OK)",
                   action="store_true", dest="strict", default=False)
parser.add_option("-n", "--blastn", help="Return BLASTN hits rather than BLASTP hits. (D: BLASTP)",
                  action="store_true", dest="blastn", default=False)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

# Do we want exact matches or LIKE %org%?
if options.strict:
    teststr = args
else:
    teststr = list('%' + s + '%' for s in args)

# We need to remove any \ in the inputs so that they dont appear in the SQL query.
# SQL doesn't understand escape characters.
teststr = [ s.replace("\\", "") for s in teststr ]

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Generate a temporary table (desiredorgs) containing the list of
# organisms for which we want to generate a CORE
#
# This query will olook like SELECT DISTINCT organism FROM processed WHERE processed.organism LIKE %org1% OR ...;
#
# ... is any number of AND statements
query = "CREATE TEMPORARY TABLE desiredorgs AS SELECT DISTINCT organism FROM processed WHERE "
for i in range(len(teststr)):
    # ! in front of an organism wildcard means NOT
    query = query + "processed.organism LIKE ? "

    if not i == len(teststr) - 1:
        query = query + "OR "
        
query = query + ";"
cur.execute(query, tuple(teststr))

# Generate a list of gene IDs to search for in the BLAST results
cur.execute("""CREATE TEMPORARY TABLE desiredgenes AS 
                 SELECT geneid FROM processed 
                 INNER JOIN desiredorgs ON desiredorgs.organism = processed.organism; """)

# Hopefully speed up the IN commands
cur.execute("""CREATE INDEX temp_desiredgenes ON desiredgenes(geneid);""")

# Generate a list of blast results with query matching one of the desiredgenes
if options.blastn:
    tbl = "blastnres_selfbit"
else:
    tbl = "blastres_selfbit"

cmd = """SELECT %s.* FROM %s
         WHERE %s.targetgene IN desiredgenes
         AND   %s.querygene  IN desiredgenes;""" %(tbl, tbl, tbl, tbl)

cur.execute(cmd)

if options.header:
    print "\t".join(header)
for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
