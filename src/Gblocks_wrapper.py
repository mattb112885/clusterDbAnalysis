#!/usr/bin/env python

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
from __future__ import print_function
import fileinput
import sys
import random
import os
import optparse
from Bio import AlignIO
from Bio import SeqIO

usage = "%prog [options] < Fasta_alignment > Fasta_alignment_filtered"
description="Run GBLOCKS with specified parameter values. GBLOCKS is a program to filter low-quality sections out of a multiple alignment"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-p", "--program", help="Name or location of your GBLOCKS program (D: Gblocks)", action="store", type="str", dest="program", default="Gblocks")
parser.add_option("-s", "--usestrict", help="Use 'strict' values according to the cited paper", action="store_true", dest="usestrict", default=False)
parser.add_option("-r", "--userelaxed", help="Use 'relaxed' values according to the cited paper (Default parameters correspond to this because the paper claims it works better for shorter alignments)", action="store_true", dest="userelaxed", default=False)
# If neither of those is specified then you will need to specify any of these that you dont want to match to the "relaxed" settings.
parser.add_option("-c", "--mincons", help="Minimum Number Of Sequences For A Conserved Position (D=9)", action="store", type="int", dest="b1", default=9)
parser.add_option("-f", "--minflank", help="Minimum Number Of Sequences For A Flank Position (D=9)", action="store", type="int", dest="b2", default=9)
parser.add_option("-n", "--maxnoncons", help="Maximum Number Of Contiguous Nonconserved Positions (D=10)", action="store", type="int", dest="b3", default=10)
parser.add_option("-m", "--minblock", help="Minimum Length Of A Block (D=5)", action="store", type="int", dest="b4", default=5)
parser.add_option("-g", "--gappos", help="Allowable gaps (D=h, up to half gaps - must be h (half), n (none), or a (any number))", action="store", type="str", dest="b5", default='h')
(options, args) = parser.parse_args()

# Sanity checks
if not options.b5 in [ 'a', 'h', 'n', 'any', 'half', 'none' ]:
    sys.stderr.write("ERROR: only valid entries for allowable gaps are a (any), n (none), or h (half)")    
if options.usestrict and options.userelaxed:
    sys.stderr.write("ERROR: Cannot use both strict and relaxed settings at the same time!\n")
    exit(2)
if options.usestrict or options.userelaxed:
    sys.stderr.write("WARNING: usestrict or userelaxed will overwrite any other parameters you have specified on the command line\n")
# For convenience
if options.b5 == 'any':
    options.b5 = 'a'
if options.b5 == 'half':
    options.b5 = 'h'
if options.b5 == 'none':
    options.b5 = 'n'

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

# Minimum Number Of Sequences For A Conserved Position
# b1 = Minimum Number Of Sequences For A Conserved Position
# b2 = Minimum Number Of Sequences For A Flank Position
# b3 = Maximum Number Of Contiguous Nonconserved Positions
# b4 = Minimum Length Of A Block
# b5 = Allowed Gap Positions
# p = y if we want a Results And Parameters File (yes, we want this)
if options.usestrict:
    b1 = 9
    b2 = 13
    b3 = 8
    b4 = 10
    b5 = 'n'
elif options.userelaxed:
    b1 = 9
    b2 = 9
    b3 = 10
    b4 = 5
    b5 = 'h'
else:
    b1 = options.b1
    b2 = options.b2
    b3 = options.b3
    b4 = options.b4
    b5 = options.b5

# Now we need to run GBlocks on this.
# Don't bother redirecting stderr but I do redirect stdout since this is a pipe command...
gblcmd = "%s %s -b1=%d -b2=%d -b3=%d -b4=%d -b5=%s -p=y > /dev/null" %(options.program, fname, b1, b2, b3, b4, b5)
sys.stderr.write("Running GBlocks with command: \n %s\n" %(gblcmd))
res = os.system(gblcmd);
# Gblocks returns 256 on success...
if not res == 256:
    sys.stderr.write("ERROR: Running GBLOCKS failed. Aborting\n")
    exit(2)

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
