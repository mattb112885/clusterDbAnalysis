#!/usr/bin/env python

# Concatenate alignments (from FASTA format)
# Alignments are all of the files that match the specified search key.
# All of the files must have exactly one representative from each organism and the
# set of organisms in each alignment must be the same.
#
# Result is exported to stdout

from Bio import SeqIO
import os
import re
import sys
import optparse

usage="%prog [options] alignment_path > concatenated_alignment"
description="""Concatinate all FASTA alignments in alignment_path.
All alignments must have exactly one represenative from each organism and all
must have representatives from the same set of organisms."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-k", "--key", help="Pattern to match (by regex) in any files you want to use, e.g. a run ID (default: use all files in alignment_path)", action="store", type="str", dest="key", default=None)
(options, args) = parser.parse_args()

if not len(args) == 1:
    sys.stderr.write("ERROR: alignment_path is a required argument (it should contain all the multiple alignments you wish to concatinate)\n")
    exit(2)

filelist = []
for filename in os.listdir(args[0]):
    # No key given --> use everything in that directory, assuming they're all what we need
    if options.key == None:
        filelist.append(filename)
    else:
        if re.search(options.key, filename) != None:
            filelist.append(filename)

if len(filelist) == 0:
    sys.stderr.write("ERROR: No files found in alignment_path or none found that match the specified key\n")
    exit(2)

# Mapping between organism ID and (aligned) sequences...
seqs = {}

for filename in filelist:
    fid = open(os.path.join(args[0], filename), "r")

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
