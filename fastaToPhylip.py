#!/usr/bin/python

# This is a pipe command.
#
# Pipe in a FASTA file
# The command will first substitute each gene ID in the fasta file with a 10-digit temporary identifier
# consistent with PHYLIP format
#

import fileinput
import sys
from Bio import AlignIO
from Bio import SeqIO

# Read the FASTA file from stdin and convert it into a phylip file
# Use list so we actually edit in-place rather than
# just editing a copy that gets destroyed later!
aln = list(AlignIO.read(sys.stdin, "fasta"))

# We will use this to convert back to the IDs in the fasta file
for i in range(len(aln)):
    newid = "S%09d" %(i)
    aln[i].id = newid

SeqIO.write(aln, sys.stdout, "phylip")
