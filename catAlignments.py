#!/usr/bin/python

# Concatenate alignments (from FASTA format)
# alignments are all of the files with extension fasta_aln
#
# Result is exported to stdout

from Bio import SeqIO
import os
import re
import sys

# Mapping between organism ID and (aligned) sequences...
seqs = {}

filelist = []
for filename in os.listdir(sys.argv[1]):
    # Search for fasta_aln files
    if re.search("\.fasta_aln$", filename) != None:
        filelist.append(filename)

for file in filelist:
    fid = open(os.path.join(sys.argv[1], filename), "r")

    records = SeqIO.parse(fid, "fasta")
    for r in records:
        myid = re.sub("\.peg\.\d+", "", r.id)
        myseq = r.seq.tostring()
        if myid in seqs:
            # concatenate sequences
            seqs[myid] = seqs[myid] + myseq
        else:
            seqs[myid] = myseq
    fid.close()

for s in seqs:
    print ">" + s
    print seqs[s]
