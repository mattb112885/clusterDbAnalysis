#!/usr/bin/python

# program [Protein_folder] [contig_folder] [out_folder] [ncores]
# Run TblastN from all the proteins in one specified FASTA file to a specified genome FASTA file
#
# I think matching this up with the clustering P/A analysis should be a separate step...

import optparse, os, re, sys
import os.path as path
from ruffus import *

usage = "%prog [Faa_folder] [Contig_folder] [Output_folder] [Options]"
description = """Run TBLASTN all vs. all for faa files in the faa_folder and fna files in the contig_folder
The program expects all protein fasta files to end with .faa or .fasta and all contig files to end with .fna """
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-t", "--threshold", help="E-value threshold (D=1E-5)", action="store", type="float", dest="threshold", default=0.00001)
parser.add_option("-n", "--ncores", help="Number of cores (D=1)", action="store", type="int", dest="ncores", default=1)
(options, args) = parser.parse_args()

if len(args) != 3:
    sys.stderr.write("ERROR: Incorrect usage \ missing arguments in tblastn_all_vs_all.py (see -h for details)\n")
    exit(2)

protDir = args[0]
genDir = args[1]
outDir = args[2]

genfilelist = []
for filename in os.listdir(path.join(os.curdir, genDir)):
    if re.search("\.fna$", filename) != None:
        genfilelist.append(filename)

protfilelist = []
for filename in os.listdir(path.join(os.curdir, protDir)):
    if re.search("\.fasta$", filename) != None or re.search("\.faa$", filename) != None:
        protfilelist.append(filename)

print genfilelist
print protfilelist

# Compile blast databases if needed
for name in genfilelist:
    try:
        handle = open(path.join(os.curdir, genDir, name) + ".nhr")
        handle.close()
    except IOError:
        cline = "makeblastdb -in %s -dbtype nucl" %(path.join(os.curdir, genDir, name))
        sys.stderr.write("%s\n" %(cline))
        os.system(cline)

params = []
for name in genfilelist:
    for name2 in protfilelist:
        outname = "%s_%s_TBLASTOUT" %(name, name2)
        try:
            handle = open(path.join(os.curdir, outDir, outname), "r")
            handle.close()
            continue
        except IOError:
            genpath = path.join(os.curdir, genDir, name)
            protpath = path.join(os.curdir, protDir, name2)
            outpath = path.join(os.curdir, outDir, outname)
            params.append([genpath, protpath, outpath, options.threshold])

if(len(params) == 0 ):
    sys.stderr.write("All requested jobs already done!\n")
    exit(0)

@parallel(params)
def runTblast(genpath, protpath, outpath, evalue):
    # Run tblastn
    # We are dealing with archaea so we want translation table 11 according to NCBI
    # Note that this will very likely choke on AUG codons translated as pyrrolysine...
    # I'm not sure exactly how to deal with those.
    #
    # Also...depending on if TBLASTN skips over STOP codons or not,
    # I might want to make the frameshift penalty very low so I can identify them.
    # However, if I have to join adjacent genes anyway for STOP codons I might as well
    # do the same thing for frameshifts...
    #
    cline = "tblastn -query %s -db %s -out %s -outfmt 6 -evalue %e -db_gencode 11" %(protpath, genpath, outpath, evalue)
    sys.stderr.write("%s\n" %(cline))
    os.system(cline)

pipeline_run([runTblast], multiprocess=options.ncores)
