#!/usr/bin/python

# THIS is a dumb genbank parser
# Biopython is too smart and doesn't let me make a fna file with multiple contigs
# if the contig name is longer than 10 characters.
#
# It is mostly intended to be used on the genbank files from RAST - it might not work
# with genbank files from other websites...
#
# Contig IDs cannot have spaces.

import fileinput
import re
import optparse

usage="%prog < genbank_file > fna file"
description="Make a contig nucleic acid FASTA file out of a genbank file. Contig names MUST have no spaces."
parser = optparse.OptionParser(usage=usage, description=description)
parser.parse_args()

spaceRemover = re.compile("\s\s+")
sequenceCleaner = re.compile("[\s\d]*")

contig = ""
seq = ""
issequence = False
for line in fileinput.input("-"):
    if line.startswith("LOCUS"):
        sub = spaceRemover.sub(" ", line)
        spl = sub.split(" ")
        contig = spl[1]
    # a line ORIGIN designates the beginning of the DNA sequence
    if line.startswith("ORIGIN"):
        issequence = True
        seq = ""
        continue
    # A line "//" designates the end of the DNA sequence...
    if line.startswith("//"):
        print "%s\t%s" %(contig, seq)
        issequence = False
        continue
    if issequence:
        # We need to remove spaces and numbers
        seq = seq + sequenceCleaner.sub("", line)
