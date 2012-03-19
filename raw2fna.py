#!/usr/bin/python

# Generates fna file from tne input file (piped) and exports it to stdout

import fileinput
for line in fileinput.input("-"):
    spl = line.strip().split('\t')
    # Protein-coding DNAs only (this also cuts off the title row)
    if not spl[2] == "peg":
        continue

    header = ">" + spl[1] + " " + spl[7]
    seq = spl[11]
    print header
    print seq
