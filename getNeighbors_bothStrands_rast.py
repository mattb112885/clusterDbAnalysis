#!/usr/bin/python

# Input: RAST raw info file (from the /raw/ folder)
#
# Outputs (up to k=12) the k-nearest neighbors of a gene on both sides. 
# Does not include circularity in general...
#
# Unlike the other version, this one will generate a neighbor table including
# BOTH strands as equivalent, rather than just ONE.

import sys
from operator import itemgetter
import re

if not len(sys.argv) == 2:
    print "Usage: ./makeNeighborTable [rast_raw_info_file]"
    exit(2)

MAXK = 5

geneTuples = []
for line in open(sys.argv[1], "r"):
    spl = line.strip().split("\t")
    # PUBSEED tables for some reason have blank lines. Lets skip those.
    if len(spl) < 6:
        continue
    # Protein-coding genes only
    if not spl[2] == "peg":
        continue

    # Gene ID contains the organism ID
    orgid = re.search("\|(\d+\.\d+)", spl[1])
    if orgid == None:
        sys.stderr.write("WARNING: No organism ID found with expected format ...\n")
        orgid = ''
    else:
        orgid = orgid.group(1)

    # Gene ID, Organism+Contig, Start, Stop, Strand
    # Add organism to contig name to avoid non-unique contig names like "contig00001" conflicting when
    # we concatinate all of them.
    geneTuples.append( (spl[1], orgid + "_" + spl[0], int(spl[4]), int(spl[5]), spl[6]) )

# Sort by contig first, then by start
sortGeneTuples = sorted(geneTuples, key=itemgetter(1,3))

# Separate them by contig
contigToTuple = {}
for tup in sortGeneTuples:
    if tup[1] in contigToTuple:
        contigToTuple[tup[1]].append(tup)
    else:
        contigToTuple[tup[1]] = [ tup ]

# Iterate through them and make a nice file
for contig in contigToTuple:
    currentTuples = contigToTuple[contig]
    for ii in range(len(currentTuples)):
        # range() function cuts off the end so we need to actually do MAXK + 1 for the rhs...
        for jj in range(-MAXK, MAXK+1):
            # No circularity - just continue if we're out of bounds.
            if ii+jj < 0:
                continue
            if ii+jj > len(currentTuples) - 1:
                continue
            # Print the gene as a neighbor - also print location and contig as a sanity check
            print "%s\t%s\t%d\t%s\t%d\t%s" %(currentTuples[ii][0], currentTuples[ii+jj][0], jj, currentTuples[ii+jj][1], currentTuples[ii+jj][2], currentTuples[ii+jj][4])
