#!/usr/bin/python
#
# BLAST_ALL_v_ALL
# Run (in parallel) blast jobs blasting all fasta files in the FASTA folder against
# all the other fasta files in that folder
#
# Fasta files must have a .fasta extension
#
# Results are outputted into FIRSTFILENAME_SECONDFILENAME.xml in the XML folder designated below.
# It will be the job of some other code to actually do the correct concatenation
#
# I may convert to table instead of XML format later...
#
# Matthew Benedict

# Run in parallel
from ruffus import *

import os
import os.path as path
import sys
import re

if len(sys.argv) != 4:
    print "Usage: Blast_all_v_all.py [Fasta_dir] [Results_dir] [Ncores]"
    exit(2)

FASTADIR = sys.argv[1]
BLASTDIR = sys.argv[2]
NCORES = sys.argv[3]

filelist = []
for filename in os.listdir(FASTADIR):
    # .fasta or .faa is acceptable
    if re.search("\.fasta$", filename) != None or re.search("\.faa$", filename) != None:
        filelist.append(filename)

for x in filelist:
    print x

# Build blast databases for all the fasta files
for target in filelist:
    try:
        handle = open(path.join(FASTADIR, target + ".phr"))
        handle.close()
    except IOError:
        cline = ["makeblastdb -in " + path.join(FASTADIR, target) + " -dbtype prot"]
        cline = "".join(cline)
        print cline
        os.system(cline)
       
params = [] 
for query in filelist:
    for target in filelist:
        # Put parallel directive here
        tmpparams = [query, target, BLASTDIR, FASTADIR]
        params.append(tmpparams)

@parallel(params)
def singleBlast(query, target, BLASTDIR, FASTADIR):
    string = query + "_" + target + "_COMBINED"
    try:
        # Don't repeat-blast; if it is already done, go to the next one.
        handle = open(path.join(BLASTDIR, string))
        handle.close()
    except IOError:
        cline = ["blastp -outfmt 6 -query " + path.join(FASTADIR, query) + " -db " + path.join(FASTADIR, target) + " -evalue 1E-5 -out " + path.join(BLASTDIR, string)]
        cline = "".join(cline)
        print cline
        os.system(cline)

pipeline_run([singleBlast], multiprocess=int(NCORES))
