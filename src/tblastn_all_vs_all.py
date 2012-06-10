#!/usr/bin/python

# program [Protein_folder] [contig_folder] [out_folder] [ncores]
# Run TblastN from all the proteins in one specified FASTA file to a specified genome FASTA file
#
# I think matching this up with the clustering P/A analysis should be a separate step...

import sys, os
import os.path as path
import re
from ruffus import *

if len(sys.argv) != 5:
    print "Usage: ./runTblastN_specGene.py [Protein_folder] [contig_folder] [output_folder] [ncores]"
    exit(2)

protDir = sys.argv[1]
genDir = sys.argv[2]
outDir = sys.argv[3]

outFmt = str(7)

genfilelist = []
for filename in os.listdir(path.join(os.curdir, genDir)):
    if re.search("\.fna$", filename) != None:
        genfilelist.append(filename)

protfilelist = []
for filename in os.listdir(path.join(os.curdir, protDir)):
    if re.search("\.fasta$", filename) != None or re.search("\.faa$", filename) != None:
        protfilelist.append(filename)


for name in genfilelist:
    try:
        handle = open(path.join(os.curdir, genDir, name) + ".nhr")
        handle.close()
    except IOError:
        cline = "makeblastdb -in " + path.join(os.curdir, genDir, name) + " -dbtype nucl"
        cline = "".join(cline)
        print cline
        os.system(cline)

params = []
for name in genfilelist:
    for name2 in protfilelist:
        outname = name + "_" + name2 + "_TBLASTOUT"
        try:
            handle = open(path.join(os.curdir, outDir, outname))
            handle.close()
            continue
        except IOError:
            genpath = path.join(os.curdir, genDir, name)
            protpath = path.join(os.curdir, protDir, name2)
            outpath = path.join(os.curdir, outDir, outname)
            params.append([genpath, protpath, outpath])

if(len(params) == 0 ):
    print "All requested jobs already done!"
    exit(0)

@parallel(params)
def runTblast(genpath, protpath, outpath):
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
    cline = "tblastn -query " + protpath + " -db " + genpath + " -out " + outpath + " -outfmt 6 -evalue 1E-5 -db_gencode 11"
    print "Running tblastn..."
    os.system(cline)

pipeline_run([runTblast], multiprocess=int(sys.argv[4]))
