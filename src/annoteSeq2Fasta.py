#!/usr/bin/python

# This is a pipe command.
# 
# Provide gene ID (column 1 or -g) annotation (column 2 or -a) and sequence (column 3 or -s) - AA or nucleotide, it doesn't matter
# I make you a FASTA file and print it to stdout. This can be used for alignments / trees / whatever.
#

import fileinput, optparse, sys

usage = "%prog [options] < gene_seq_table > fasta_file"
description = "Turn a table containing gene IDs, annotations, and sequences into a FASTA file"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--geneid", help="Column number (start from 1) for gene id", action="store", type="int", dest="geneidcol", default=1)
parser.add_option("-a", "--annote", help="Column number (start from 1) for annotation / FASTA header (optional, by default no annotation is included)", action="store", type="int", dest="annotecol", default=None)
parser.add_option("-s", "--seqcol", help="Column number (start from 1) for sequence column", action="store", type="int", dest="seqcol", default=3)
(options, args) = parser.parse_args()

geneidcol = options.geneidcol - 1
seqcol = options.seqcol - 1
if options.annotecol is not None:
    annotecol = options.annotecol - 1
else:
    annotecol = None

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    geneid = spl[geneidcol]
    seq = spl[seqcol]
    if annotecol is None:
        annote = ""
    else:
        annote = " %s" %(spl[annotecol])
    print ">%s%s\n%s" %(geneid, annote, seq)
