#!/usr/bin/python

# Returns a list of all bidirectional best hits 
# between any pair of organisms from the
# blast results table.

# The best hit is defined by E-value but could later be defined by any
# other score we want with appropriate input arguments...

import sqlite3
import optparse
import sys
import math

usage="%prog [options] > BBH_table"
description = "Return a list of bidirectional best blast hits based on the specified scoring criteria. Output table has (tab-delimited): Query gene, target gene, query genome, forward score, backward score"
parser = optparse.OptionParser(description=description, usage=usage)
parser.add_option("-m", "--method", help="Scoring metric to use to define best hit (D=evalue)", action="store", type="str", dest="method", default="evalue")
(options, args) = parser.parse_args()

okmethods = [ "evalue", "maxbit", "minbit" ]

if not options.method in okmethods:
    print "ERROR: Specified method in db_bidirectionalBestHits.py not implemented."
    print "Implemented methods: %s" %("\t".join(okmethods))
    exit(2)

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Get a list of BLAST results with
# organism for query and target attached.
query = """ SELECT blastres_selfbit.*, p1.organism as queryorg, p2.organism AS targetorg
            FROM blastres_selfbit
            INNER JOIN processed as p1 ON p1.geneid = blastres_selfbit.querygene
            INNER JOIN processed as p2 ON p2.geneid = blastres_selfbit.targetgene; """

sys.stderr.write("Executing SQL query...\n")
cur.execute(query)

# Link from (query gene, target organism) to (target gene, score, query organism)
Best_pairs = {}
sys.stderr.write("Reading best hits and calculating scores... this could take some time\n")
n = 0
for s in cur:

    if n - (n/100000)*100000 == 0:
        sys.stderr.write("%d\n" %(n) )

    ls = [ str(k) for k in s ]
    # (Query gene, target organism)
    mypair = ( ls[0], ls[15] )
    # Smaller E-values are better (hence the - sign)
    if options.method == "evalue":
        myscore = -math.log10(float(ls[10]) + 1E-200)
    # Bigger maxbit and minbit are better
    elif options.method == "maxbit":
        myscore = float(ls[11])/max(float(ls[12]), float(ls[13]))
    elif options.method == "minbit": 
        myscore = float(ls[11])/min(float(ls[12]), float(ls[13]))

    if mypair in Best_pairs:
        pairpair = Best_pairs[mypair]
        if pairpair[1] < myscore:
            Best_pairs[mypair] = (ls[1], myscore, ls[14])
    else:
        Best_pairs[mypair] = (ls[1], myscore, ls[14])
    n += 1

# Look at all the best queries and see if the target gene paired with the same organism
# has a best hit as the query.
sys.stderr.write("Calculating bidirectional bests...\n")
for pair in Best_pairs:
    forward_hit = Best_pairs[pair]

    # We want the highest hit TO THE QUERY GENOME.
    back_pair = ( forward_hit[0], forward_hit[2] )

    if not back_pair in Best_pairs:
        sys.stderr.write("WARNING: Hit %s,%s,%s,%1.4f had no backwards hit\n" %(pair[0], forward_hit[0], pair[1], forward_hit[1]))
        continue
    backward_hit = Best_pairs[ back_pair ]
    # Bidirectional best: The best hit in the forward and backward direction was the same!
    # I print the forward and backward scores because E-values aren't necessarily symmetric...
    if backward_hit[0] == pair[0]:
        # Query gene, target gene, query genome, forward score, backward score
        print "%s\t%s\t%s\t%1.4f\t%1.4f" %(pair[0], forward_hit[0], pair[1], forward_hit[1], backward_hit[1])

con.close()
