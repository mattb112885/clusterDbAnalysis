#!/usr/bin/python

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

if not os.path.exists(target_db + ".phr"):
    # The -S 1 is because NCBI recommended NOT using scaling even
    # though they set it as the default anyway.
    # See:
    #  ftp://ftp.ncbi.nih.gov/blast/documents/rpsblast.html
    # This is the docs for the old version - they removed the warning
    # from the new version but didn't explain why...
    os.system("formatrpsdb -i \"%s\" -S 1" %(target_db))

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
    cmd = "rpsblast -i %s -d %s -o %s -m 8 -e %e -P 1;" %(faa, target_db, outfile, evalue)
    sys.stderr.write("%s\n" %(cmd))
    os.system(cmd)

pipeline_run([singleRpsBlast], multiprocess=options.numcores)
