#!/usr/bin/python
#
# BLAST_ALL_v_ALL
# Run (in parallel) blast jobs blasting all fasta files in the FASTA folder against
# all the other fasta files in that folder
#
# Fasta files must have a .fasta extension or a .faa extension (to prevent trying to blast against the database files)
#
# Results are outputted in a standard name combining the names of the two fasta files used
# to generate them. Format is -m9 (or -outfmt 6 with the new BLAST) - tab-delimited with default fields
# and no comments.
#
# Matthew Benedict

# Run in parallel
from ruffus import *

import os
import os.path as path
import sys
import re
import optparse

usage="%prog [-n] Fasta_dir Results_dir Ncores"
description="""Automatically make blast databases for all files in fasta_dir 
(must have .faa or .fasta extension) in parallel using ncores processers, and 
store results in results_dir"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--blastn", help="Perform BLASTN rather than BLASTP (D: BlastP). Runs actual BLASTN, not MEGABLAST.", 
                  action="store_true", dest="blastn", default=False)
(options, args) = parser.parse_args()

if len(args) != 3:
    sys.stderr.write("ERROR: Fasta_dir, results_dir and ncores are all required. Type Blast_all_vs_all.py -h for details\n")
    exit(2)

FASTADIR = args[0]
BLASTDIR = args[1]
NCORES = args[2]

filelist = []
for filename in os.listdir(FASTADIR):
    # .fasta or .faa is acceptable
    if re.search("\.fasta$", filename) != None or re.search("\.faa$", filename) != None or re.search("\.fna$", filename) != None:
        filelist.append(filename)

for x in filelist:
    sys.stderr.write("%s\n" %(x))

if len(filelist) == 0:
    sys.stderr.write("""ERROR: No files with .faa or .fasta extension found in specified fasta directory %s. Aborting..\n"""
                     %(FASTADIR) )
    exit(2)

# Build blast databases for all the fasta files
for target in filelist:
    try:
        if options.blastn:
            handle = open(path.join(FASTADIR, target + ".nhr"))
            handle.close()
        else:
            handle = open(path.join(FASTADIR, target + ".phr"))
            handle.close()
    except IOError:
        if options.blastn:
            cline = "makeblastdb -in %s -dbtype nucl" %(path.join(FASTADIR, target))
        else:
            cline = "makeblastdb -in %s -dbtype prot" %(path.join(FASTADIR, target))
        sys.stderr.write("%s\n" %(cline))
        os.system(cline)
       
params = [] 
for query in filelist:
    for target in filelist:
        # Put parallel directive here
        tmpparams = [query, target, BLASTDIR, FASTADIR, options]
        params.append(tmpparams)

@parallel(params)
def singleBlast(query, target, BLASTDIR, FASTADIR, options):
    string = "%s_%s_COMBINED" %(query, target)
    try:
        # Don't repeat-blast; if it is already done, go to the next one.
        handle = open(path.join(BLASTDIR, string))
        handle.close()
    except IOError:
        if options.blastn:
            # I use a higher cutoff here since it doesn't take much to store \ query blastn results anyway.
            cline = "blastn -task blastn -outfmt 6 -query %s -db %s -evalue 1 -out %s" %( path.join(FASTADIR, query), path.join(FASTADIR, target), path.join(BLASTDIR, string) )
        else:
            cline = "blastp -outfmt 6 -query %s -db %s -evalue 1E-5 -out %s" %( path.join(FASTADIR, query), path.join(FASTADIR, target), path.join(BLASTDIR, string) )
        sys.stderr.write("%s\n" %(cline))
        os.system(cline)

pipeline_run([singleBlast], multiprocess=int(NCORES))
