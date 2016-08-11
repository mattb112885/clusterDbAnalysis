#!/usr/bin/env python

# Generates fna file from tne input file (piped) and exports it to stdout

import fileinput
import optparse

usage = "%prog < RAW_file > faa_file"
description="Make an amino-acid fasta file for all the proteins in the piped-in RAW file"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split('\t')
    # Protein-coding DNAs only (this also cuts off the title row)
    if not spl[2] == "peg":
        continue

    header = ">" + spl[1] + " " + spl[7]
    seq = spl[12]
    print(header)
    print(seq)
