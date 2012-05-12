#!/usr/bin/python

# This is a pipe command.
#
# Pipe in a FASTA file
# The command will first substitute each gene ID in the fasta file with a 10-digit temporary identifier
# consistent with PHYLIP format
#
# Then the file converts the fasta file to a temporary PHYLIP file.
# After running RAXML with that temporary PHYLIP file,
# the program will take the resulting .nwk file and convert the IDs back.
#
# I wrote this because the svr_tree function didn't work for me when trying to run it
# with RAXML, and I need something...

import fileinput
import sys
import random
import os
from Bio import AlignIO
from Bio import SeqIO
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-b", "--bootstraps", help="Number of bootstraps (D=0)", action="store", type="int", dest="NUMBOOTS", default=0)
parser.add_option("-T", "--numthreads", help="Number of threads, must be more than 1 (D=2)", action="store", type="int", dest="NTHREADS", default=2)
parser.add_option("-k", "--nocleanup", help="Set this flag to keep intermediate RAXML and PHYLIP files (D=false, delete these files)", action="store_false", dest="CLEANUP", default=True)
parser.add_option("-m", "--model", help="Specify model to use with RAXML (D=PROTGAMMAWAG)", action="store", type="str", dest="MODEL", default="PROTGAMMAWAG")
parser.add_option("-p", "--program", help="Specify the name of the RAXML program to use (D=raxmlHPC-PTHREADS)", action="store", type="str", dest="PROGRAM", default="raxmlHPC-PTHREADS")

(options, args) = parser.parse_args()

NUMBOOTS=options.NUMBOOTS
NUMTHREADS=options.NTHREADS
CLEANUP=options.CLEANUP
MODEL=options.MODEL
PROGRAM=options.PROGRAM

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

# Make a temporary random file to place the phylip (must have write permission in the current directory
# I would've used /tmp/ but RAXML assumes you want the files in the current directory and
# gets very angry when you end up specifying [currentdirectory]/./tmp/...
#
# For whatever reason, there are no problems with the phylip writing in SeqIO
# but it doesn't work in AlignIO, while for reading the FASTA it's the opposite case.
# Whatever.
#
# This isn't exactly thread-safe
fname = "%d.phi" %(random.randint(0,2**30))

fid = open(fname, "w")
SeqIO.write(aln, fid, "phylip")
fid.close()

# Now we need to run RAXML on this.
# I should set up some options for people to pass in here
# But for now I'll just do this
INFILE=fname
OUTFILE="%s.nwk" %(fname)

# Required part
raxcmd = "%s -s %s -n %s -m %s -T %d -p 123456" %(PROGRAM, INFILE, OUTFILE, MODEL, NUMTHREADS)

# Optional part
if NUMBOOTS > 0:
    raxcmd = "%s -f a -x 123456 -# %d" %(raxcmd, NUMBOOTS)

# Since this is a pipe command and all this stuff is saved in one of the output args anyway
# we don't want the outputs if they are sent to stdin (stderr is fine).
raxcmd = "%s > /dev/null" %(raxcmd)

# Run RAXML
sys.stderr.write("Running RaxML with command: \n %s\n" %(raxcmd))
os.system(raxcmd);

# Get the best tree from those results and replace the names
# with the real names
treestr = "".join([ line.strip() for line in open("RAxML_bestTree.%s" %(OUTFILE)) ])

for sub in subToReal:
    treestr = treestr.replace(sub, subToReal[sub])

print treestr

if CLEANUP:
    # * to remove the "reduced" files if any are created by RAXML
    os.system("rm %s*" %(INFILE) )
    os.system("rm RAxML_*.%s" %(OUTFILE))

