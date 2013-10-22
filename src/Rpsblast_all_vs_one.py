#!/usr/bin/env python

import optparse, os, sys
from ruffus import *

usage = "%prog [options] target_db inputdir outputdir"
description = """Parallelize RPSBLAST by running each organism in the inputdir
separately against the same target database."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--numcores", help="Number of cores to use (D: 1)", 
                  action="store", type="int", dest="numcores", default=1)
parser.add_option("-e", "--extension", help="Extension for input fasta files (D: faa)", 
                  action="store", type="str", dest="extension", default="faa")
parser.add_option("-c", "--cutoff", help="E-value cutoff (D: 1E-5)", 
                  action="store", type="float", dest="cutoff", default=1E-5)
(options, args) = parser.parse_args()

if len(args) != 3:
    sys.stderr.write("ERROR: target_db, inputdir, and outputdir are necessary input arguments\n")
    exit(2)

target_db = args[0]
inputdir = args[1]
outputdir = args[2]

faa_list = []
for f in os.listdir(os.path.join(args[1])):
    if f.endswith(".%s" %(options.extension)):
        faa_list.append( os.path.join(args[1], f) )

evalue = options.cutoff

arglist = []
for faa in faa_list:
    arglist.append( [ faa, target_db, outputdir, evalue ] )

@parallel(arglist)
def singleRpsBlast(faa, target_db, outputdir, evalue):
    outfile = os.path.join(outputdir, "%s_rpsout" %(os.path.basename(faa)))
    if os.path.exists(outfile):
        return
    cmd = "rpsblast -query %s -db %s -out %s -outfmt 6 -evalue %e" %(faa, target_db, outfile, evalue)
    # This was the command we ran on the old BLAST
    # NCBI updated the RPSBLAST with recent BLAST+ and the CDD no longer works with the old compiler
    # so... lets just make everything conform to thh new standard.
#    cmd = "rpsblast -i %s -d %s -o %s -m 8 -e %e -P 1;" %(faa, target_db, outfile, evalue)
    sys.stderr.write("%s\n" %(cmd))
    os.system(cmd)

pipeline_run([singleRpsBlast], multiprocess=options.numcores)
