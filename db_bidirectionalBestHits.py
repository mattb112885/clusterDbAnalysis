#!/usr/bin/python

# Returns a list of all bidirectional best hits 
# between any pair of organisms from the
# blast results table.

# The best hit is defined by E-value but could later be defined by any
# other score we want with appropriate input arguments...

import sqlite3
import optparse
import sys

description = "Return a list of bidirectional best blast hits based on the specified scoring criteria"
parser = optparse.OptionParser(description=description)
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
# organism for target attached.
# We don't need to bother with the query's organism...
query = """ SELECT blastres_selfbit.*, processed.organism AS targetorg
            FROM blastres_selfbit
            INNER JOIN processed ON geneid = blastres_selfbit.targetgene; """

sys.stderr.write("Executing SQL query...\n")
cur.execute(query)

# Link from (query gene, target organism) to (target gene, score)
Best_pairs = {}
sys.stderr.write("Reading best hits and calculating scores... this could take some time\n")
n = 0
for s in cur:

    if n - (n/10000)*10000 == 0:
        print n

    ls = [ str(k) for k in s ]
    # (Query gene, target organism)
    mypair = ( ls[0], ls[14] )

    # Smaller E-values are better (hence the - sign)
    if options.method == "evalue":
        myscore = -float(ls[10])
    # Bigger maxbit and minbit are better
    elif options.method == "maxbit":
        myscore = float(ls[11])/max(float(ls[12]), float(ls[13]))
    elif options.method == "minbit": 
        myscore = float(ls[11])/min(float(ls[12]), float(ls[13]))

    if mypair in Best_pairs:
        pairpair = Best_pairs[mypair]
        if pairpair[1] < myscore:
            Best_pairs[mypair] = (ls[1], myscore)
    else:
        Best_pairs[mypair] = (ls[1], myscore)
    n += 1

# Look at all the best queries and see if the target gene paired with the same organism
# has a best hit as the query.
sys.stderr.write("Calculating bidirectional bests...\n")
for pair in Best_pairs:
    forward_hit = Best_pairs[pair]
    backward_hit = Best_pairs(forward_hit[0], pair[1])
    # Bidirectional best: The best hit in the forward and backward direction was the same!
    # I print the forward and backward scores because E-values aren't necessarily symmetric...
    if backward_hit[0] == pair[0]:
        print "%s\t%s\t%s\t%1.4f\t%1.4f" %(pair[0], forward_hit[0], pair[1], forward_hit[1], backward_hit[1])

con.close()
