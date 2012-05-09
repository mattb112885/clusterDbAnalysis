#!/usr/bin/python

# Deals with more too-long identifier name bullshit from other programs
#
# Provide it both an input and an output fasta file name and it will go ahead
# and shorten the names for you and call GBlocks (with "relaxed" parameters that are better
# for more closely-related or shorter alignments) to trim the alignment to high-quality sites.
#
# Default parameters are described in the paper
# "Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments"
# C.F.:
# http://sysbio.oxfordjournals.org/content/56/4/564.full
#
import sys
import random
import os
from Bio import AlignIO
from Bio import SeqIO

if not len(sys.argv) == 3:
    sys.stderr.write("Usage: ./Gblocks_wrapper.py [input alignment] [output alignment]\n")
    exit(2)

# Read the FASTA file from stdin and convert it into a new fasta file
aln = list(AlignIO.read(open(sys.argv[1], "r"), "fasta"))

# We will use this to convert back to the IDs in the fasta file
subToReal = {}
for i in range(len(aln)):
    newid = "S%09d" %(i)
    subToReal[newid] = aln[i].description
    aln[i].id = newid
    # Gblocks STILL complains that the name is too fucking long even if the ID is shortened. FINE, I'll throw out the fucking annotation.
    aln[i].description = ""

fname = "%d.faa" %(random.randint(0,2**30))

fid = open(fname, "w")
SeqIO.write(aln, fid, "fasta")
fid.close()

# Now we need to run GBlocks on this.
gblcmd = "Gblocks %s -b1=9 -b2=9 -b3=10 -b4=5 -b5=h -p=y > /dev/null 2> /dev/null" %(fname)
sys.stderr.write("Running GBlocks with command: \n %s\n" %(gblcmd))
os.system(gblcmd);

# Open up the -gb file created by Gblocks
# and put the real IDs back in there.
# SIGH...
aln = list(AlignIO.read(open(fname + "-gb", "r"), "fasta"))
for i in range(len(aln)):
    des = subToReal[aln[i].id]
    aln[i].id=des
    aln[i].description=""

fid = open(sys.argv[2], "w")
SeqIO.write(aln, fid, "fasta")
fid.close()

# Clean up temporary files
rmcmd = "rm %s*" %(fname)
os.system(rmcmd)
