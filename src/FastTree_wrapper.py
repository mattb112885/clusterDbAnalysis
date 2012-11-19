#!/usr/bin/python

# This is a pipe command.
#
# Pipe in a FASTA file
#
# The command will first substitute each gene ID in the fasta file with a 10-digit temporary identifier
# consistent with PHYLIP format (needed for global bootstrapping)
#
# Then the file converts the fasta file to a temporary PHYLIP file.
# After running RAXML with that temporary PHYLIP file,
# the program will take the resulting .nwk file and convert the IDs back.
#
# This is meant to abstract away some of the details of using FASTTREE for 
# choosing a model, and especially for performing global boostrap analysis
# using FASTTREE

import fileinput
import sys
import random
import os
from Bio import AlignIO
from Bio import SeqIO
from optparse import OptionParser

description = "Wrapper for FASTTREE to allow for global bootstrapping in addition to normal functionality"
usage="%prog [options] < FASTA_file > Newick_file"
parser = OptionParser(description=description, usage=usage)
parser.add_option("-b", "--bootstraps", help="Number of bootstraps (D=none)", action="store", type="int", dest="NUMBOOTS", default=0)
parser.add_option("-g", "--globalboots", help="Perform global bootstrapping analysis. Requires Phylip with SEQBOOT and CompareToBootstrap.pl from the FASTTREE doc page (D:Local bootstrap only)", action="store_true", dest="globalboots", default=False)
parser.add_option("-n", "--nogamma", help="Do not apply gamma20 likelihood calculation", action="store_true", dest="nogamma", default=False)
parser.add_option("-m", "--model", help="Specify model to use with FASTTREE (D=WAG)", action="store", type="str", dest="MODEL", default="wag")
parser.add_option("-p", "--program", help="Specify the name of the FASTTREE program to use - it must be in your PATH (D=FastTreeMP)", action="store", type="str", dest="PROGRAM", default="FastTreeMP")
parser.add_option("-t", "--nucleotides", help="Set this flag if you are inputting a nucleotide alignment (D = Protein alignment)", action="store_true", dest="nt", default=False)

(options, args) = parser.parse_args()

##
# List of models used in FastTree
# Note jtt is the default for AAs and jc for NTs so if we get that we just won't add a model flag
##
okmodels = [ "wag", "jtt", "jc", "gtr"]

NUMBOOTS=options.NUMBOOTS
GLOBALBOOTS=options.globalboots
MODEL=options.MODEL.lower()
PROGRAM=options.PROGRAM
NUCLEOTIDE=options.nt
USEGAMMA=not options.nogamma

if MODEL not in okmodels:
    sys.stderr.write("ERROR: specified model not supported by FastTree\n")
    exit(2)

if GLOBALBOOTS and NUMBOOTS == 0:
    sys.stderr.write("ERROR: Specifying global bootstrap without any bootstraps doesnt make sense! Did you forget to specify -b?\n")
    exit(2)

# Read the FASTA file from stdin and convert it into a phylip file
# Use list so we actually edit in-place rather than
# just editing a copy that gets destroyed later!
aln = list(AlignIO.read(sys.stdin, "fasta"))

# We will use this to convert back to the IDs in the fasta file
subToReal = {}
for i in range(len(aln)):
    newid = "S%09d" %(i)
    subToReal[newid] = aln[i].id
    aln[i].id = newid

#############
# Make a temporary random file to place the phylip (must have write permission in the current directory)
#
# For whatever reason, there are no problems with the phylip writing in SeqIO
# but it doesn't work in AlignIO, while for reading the FASTA it's the opposite case.
# Whatever.
#
# This isn't exactly thread-safe
##############
fname = "%d.phi" %(random.randint(0,2**30))

fid = open(fname, "w")
SeqIO.write(aln, fid, "phylip")
fid.close()

##############
# Set up first FASTTREE run (and only one without global bootstraps)
##############

# Specify nucleotides
ntflag = ""
if NUCLEOTIDE:
    # I let FASTTREE tell the user if the model doesn't work with the specified character type.
    ntflag = "-nt"

# Flag to set for model
modelflag = ""
if MODEL=="wag":
    modelflag = "-wag"
elif MODEL=="gtr":
    modelflag = "-gtr"
elif MODEL=="jtt" or MODEL=="jc":
    # Defaults
    pass;

# Specify number of bootstrap values
#
# Don't bother with this for the global since we'll be replacing the values anyway
# with something from the global bootstrapping.
bootflag = ""
if GLOBALBOOTS:
    bootflag = "-boot 0"
else:
    bootflag = "-boot %d" %(NUMBOOTS)

# Need the gamma log file if CONSEL is to be used; use gamma anyway if people want it.
gammaflag = ""
if USEGAMMA:
    gammaflag = "-gamma"

# Initial command to run FASTTREE
FastTreeCmd = "cat %s | %s %s %s %s %s > %s.nwk" %(fname, PROGRAM, ntflag, modelflag, bootflag, gammaflag, fname)
sys.stderr.write(FastTreeCmd + "\n")
os.system(FastTreeCmd)

# Deal with global bootstrap if requested
if GLOBALBOOTS:
    # Generate the phylip bootstrap file
    os.system("phylipSeqbootScript.sh %s %s_seqboot %d" %(fname, fname, NUMBOOTS))
    # Run FastTree with the same parameters but on the bootstrap alignments
    FastTreeCmd = "cat %s_seqboot | %s %s %s %s %s -n %d > %s_seqboot.nwk" %(fname, PROGRAM, ntflag, modelflag, bootflag, gammaflag, NUMBOOTS, fname)
    sys.stderr.write(FastTreeCmd + "\n")
    os.system(FastTreeCmd)
    # Run the CompareToBootstrap.pl program needed to do something with the global bootstrap values...
    os.system("CompareToBootstrap.pl -tree %s -boot %s_seqboot.nwk > %s_global.nwk" %(fname, fname, fname))
    os.system("mv %s_global.nwk %s.nwk" %(fname, fname) )

treestr = "".join([ line.strip('\r\n') for line in open("%s.nwk" %(fname)) ])
# Substitute names back for the final tree output
for sub in subToReal:
    treestr = treestr.replace(sub, subToReal[sub])

# Print the final tree to stdout
print treestr

# Clean up temporary files
#os.system("rm %s*" %(fname) )
