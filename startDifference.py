#!/usr/bin/python

# Calcualtes the minimum distance (negatives indicate overlap) between
# consecutive genes
#
# This is a pipe function. Pipe in the RAW file

import fileinput
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
    spl = line.strip().split("\t")
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
    print str(minval) + "\t" + lastgene + "\t" + lastannote + "\t" +  gene + "\t" + annote

    laststart = int(spl[4])
    laststop = int(spl[5])
    lastgene = spl[1]
    lastannote = spl[7]
