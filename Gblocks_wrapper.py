#!/usr/bin/python

# This is a pipe command.
#
# Pipe in an alignment
# (e.g. from MAFFT, clustal, whatever)
# in FASTA format
#
# Outputs the results of Gblocks to stdout.
#
# Default parameters are "relaxed" parameters described in the paper
# "Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments"
# C.F.:
# http://sysbio.oxfordjournals.org/content/56/4/564.full
#
# According to the paper the relaxed parameters are better when dealing with short alignments (i.e. only one gene)
#
import fileinput
import sys
import random
import os
from Bio import AlignIO
from Bio import SeqIO

# Read the FASTA file from stdin and convert it into a new fasta file
aln = list(AlignIO.read(sys.stdin, "fasta"))

# We will use this to convert back to the real IDs later
subToReal = {}
for i in range(len(aln)):
    newid = "S%09d" %(i)
    subToReal[newid] = aln[i].description
    aln[i].id = newid
    # Gblocks STILL complains that the name is too long even if the ID is shortened. FINE, I'll throw out the annotation too...
    aln[i].description = ""

fname = "%d.faa" %(random.randint(0,2**30))

fid = open(fname, "w")
SeqIO.write(aln, fid, "fasta")
fid.close()

# Now we need to run GBlocks on this.
# Don't bother redirecting stderr but I do redirect stdout since this is a pipe command...
gblcmd = "Gblocks %s -b1=9 -b2=9 -b3=10 -b4=5 -b5=h -p=y > /dev/null" %(fname)
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

SeqIO.write(aln, sys.stdout, "fasta")

# Clean up temporary files
rmcmd = "rm %s*" %(fname)
os.system(rmcmd)
