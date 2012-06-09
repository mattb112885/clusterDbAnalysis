#!/usr/bin/python

# Calcualtes the minimum distance (negatives indicate overlap) between
# consecutive genes
#
# This is a pipe function. Pipe in the RAW file

#
# FIXME : If I want to keep this function around I should have it
# sort things and identify if they're in the same contig and such...
import fileinput
import optparse

usage = "%prog < RAW_file > Min_distance_file"
description="Calculates the minimum distance between consecutive genes. Negative values indicaite overlaps. Does not look at whether contigs are the same or not."
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

firstone = True

start = None
stop = None
gene = None
annote = None
laststart = None
laststop = None
lastgene = None
lastannote = None

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    if firstone:
        laststart = int(spl[4])
        laststop = int(spl[5])
        lastgene = spl[1]
        lastannote = spl[7]
        firstone = False
        continue

    start = int(spl[4])
    stop = int(spl[5])
    gene = spl[1]
    annote = spl[7]

    minval = min( start - laststart, start - laststop, stop - laststart, stop - laststop)
    print "%d\t%s\t%s\t%s\t%s" %(minval, lastgene, lastannote, gene, annote)

    laststart = int(spl[4])
    laststop = int(spl[5])
    lastgene = spl[1]
    lastannote = spl[7]
