#!/usr/bin/python

# Concatenate alignments (from FASTA format)
# alignments are all of the files with extension fasta_aln
#
# Result is exported to stdout

from Bio import SeqIO
import os
import re
import sys

if not len(sys.argv) == 2:
    print "Usage: ./catAlignments.py [alignment path]"
    exit(2)

filelist = []
for filename in os.listdir(sys.argv[1]):
    # Search for fasta_aln files
    if re.search("\.fasta_aln$", filename) != None:
        filelist.append(filename)

if len(filelist) == 0:
    print "ERROR: No files found with expected extension fasta_aln"
    exit(2)

# Mapping between organism ID and (aligned) sequences...
seqs = {}

for filename in filelist:
    fid = open(os.path.join(sys.argv[1], filename), "r")

    records = SeqIO.parse(fid, "fasta")
    for r in records:
        myid = re.sub("\.peg\.\d+", "", r.id)
        if myid in seqs:
            # Add on to existing sequence
            seqs[myid] = seqs[myid] + str(r.seq)
        else:
            seqs[myid] = str(r.seq)

    fid.close()

# Print concatinated result
for s in seqs:
    print ">" + s
    print seqs[s]
