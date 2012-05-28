#!/usr/bin/python

# This is a pipe command.
# 
# Provide gene ID (column 1 or -g) annotation (column 2 or -a) and sequence (column 3 or -s) - AA or nucleotide, it doesn't matter
# I make you a FASTA file and print it to stdout. This can be used for alignments / trees / whatever.
#

import fileinput, optparse, sys

parser = optparse.OptionParser()
parser.add_option("-g", "--geneid", help="Column number (start from 1) for gene id", action="store", type="int", dest="geneidcol", default=1)
parser.add_option("-a", "--annote", help="Column number (start from 1) for annotation / FASTA header", action="store", type="int", dest="annotecol", default=2)
parser.add_option("-s", "--seqcol", help="Column number (start from 1) for sequence column", action="store", type="int", dest="seqcol", default=3)
(options, args) = parser.parse_args()

geneidcol = options.geneidcol - 1
annotecol = options.annotecol - 1
seqcol = options.seqcol - 1

for line in fileinput.input("-"):
    spl = line.strip('\n').split("\t")
    geneid = spl[geneidcol]
    annote = spl[annotecol]
    seq = spl[seqcol]

#    if len(geneid) > 10:
#        sys.stderr.write("WARNING: Gene ID lengths greater than 10 identified. This could cause problems for PHYLIP format-dependent programs like RAXML\n")

    print ">" + geneid + " " + annote
    print seq
