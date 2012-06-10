#!/usr/bin/python

# Run from root directory
# CONCLUSION: The only genes RAST returns with internal STOP
# codons are the pyrrolysine-containing genes (i.e. it gets that correct).
# Thus we can check for the internal STOPs
# [seen here]
# and we can see if any additional STOPs aside from those present in these
# genes are added to the TBLASTN hits

import sqlite3, optparse

usage = "%prog > [genes_with_internal_stops]"
description="Find genes with internal stops that are present in the database. Internal stops are defined as internal TAG, TGA, or TAA codons"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# STOP codons
# Some of these can be skipped over
# What I want to do is identify which of the SEED
# genes (if any) have internal in-frame stop codons.
STOP = ["tga", "tag", "taa"]

cur.execute('select geneid, nucseq, annotation from rawdata where ftype = "peg";')

for rec in cur:
    l = list(rec)
    myid = l[0]
    myseq = l[1]
    annote = l[2]

    currentIdx = 0
    # Skip over the last codon because that's obviously going to be a STOP
    while currentIdx < len(myseq) - 4:
        currentCodon = myseq[currentIdx:currentIdx+3]
        if currentCodon in STOP:
            print myid + "\t" + currentCodon + "\t" + annote
        currentIdx = currentIdx + 3


con.close()
