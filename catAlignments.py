#!/usr/bin/python

# Concatenate alignments (from FASTA format)
# alignments are all of the files with extension fasta_aln
#
# Result is exported to stdout

from Bio import SeqIO
import os
import re
import sys

if not len(sys.argv) == 3:
    print "Usage: ./catAlignments.py [alignment path] [Searchkey]"
    print "Searchkey is a pattern in the files to search for."
    print "It is generally the first part of a run ID for which you want to generate"
    print "a tree, e.g. 8XCM"
    exit(2)

filelist = []
for filename in os.listdir(sys.argv[1]):
    # Search for files with the keyword in them and with at least fasta_aln (maybe fasta_aln_trimmed)
    if re.search("fasta_aln") != None and re.search(sys.argv[2]) != None:
        filelist.append(filename)

if len(filelist) == 0:
    print "ERROR: No files found with expected extension fasta_aln and search key Searchkey"
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
