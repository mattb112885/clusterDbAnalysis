#!/usr/bin/env python

# Returns a list of all bidirectional best hits 
# between any pair of organisms from the
# blast results table.

# The best hit is defined by E-value but could later be defined by any
# other score we want with appropriate input arguments...

import fileinput
import math
import optparse
import sqlite3
import sys
from FileLocator import *
from ClusterFuncs import *

okmethods = [ "evalue", "maxbit", "minbit" ]

usage="""%prog [options] > BBH_table"""
description = "Return a list of bidirectional best blast hits based on the specified scoring criteria. Output table has (tab-delimited): Query gene, target gene, query genome, forward score, backward score"
parser = optparse.OptionParser(description=description, usage=usage)
parser.add_option("-m", "--method", help="Scoring metric to use to define best hit (D=evalue). Defined methods: %s" %(" ".join(okmethods)),
                  action="store", type="str", dest="method", default="evalue")
parser.add_option("-r", "--runid", help="Get bidirectional best BLAST hits for organisms in this cluster run only (D: Get them for all organisms in the database)",
                  action="store", type="str", dest="runid", default=None)
parser.add_option("-f", "--orgfile", help="File containing s list of organisms to which to limit search. Use \"-\" for stdin. Cannot use with -r. (D: Get hits between all organisms in the database)",
                  action="store", type="str", dest="orgfile", default=None)
parser.add_option("-o", "--orgcol", help="Column number for organism ids (ignored unless -f is specified), starting from 1 (D=1)", action="store", type="int", dest="oc", default=1)
(options, args) = parser.parse_args()

oc = options.oc - 1

if not options.method in okmethods:
    sys.stderr.write("ERROR: Specified method in db_bidirectionalBestHits.py not implemented.\n")
    sys.stderr.write("Implemented methods: %s\n" %("\t".join(okmethods)))
    exit(2)

if options.runid is not None and options.orgfile is not None:
    sys.stderr.write("ERROR: runid (-r) and orgfile (-f) flags are incompatible - pick one or the other.\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

orglist = None
if options.runid is not None:
    orglist = set(getOrganismsInClusterRun(options.runid, cur))
elif options.orgfile is not None:
    orglist = set()
    for line in fileinput.input(options.orgfile):
        spl = line.strip("\r\n").split("\t")
        orglist.add(spl[oc])

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
    # Filter for specific organism pairings
    if orglist is not None:
        if str(s[14]) not in orglist or str(s[15]) not in orglist:
            continue

    # This is just a counter.
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

if len(list(Best_pairs.keys())) == 0:
    sys.stderr.write("ERROR: No bidirectional best hits found (were the organisms passed into this function correct?)\n")
    exit(2)

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
        print("%s\t%s\t%s\t%1.4f\t%1.4f" %(pair[0], forward_hit[0], pair[1], forward_hit[1], backward_hit[1]))

con.close()
