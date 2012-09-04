#!/usr/bin/python

# Wrapper script for TBLASTN verification of presence\absence
# across different genomes

import fileinput, os, optparse, random, sqlite3, sys
from locateDatabase import *

usage = "%prog (-d|-f|-o) contig_inputs [options] < Protein_ids > Tblastn_table"
description = """Attempts to run TBLASTN and identify
missing genes. It identifies called genes that match the hit location and
also tries to find genes on the opposite strand that conflict.

Either you must specify the compiled contig database or must have run ./main3.sh and
specify a single organism or list of organism IDs against which to perform the BLAST."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-o", "--organism", help="Organism ID to BLAST against.", action="store", type="str", dest="org", default=None)
parser.add_option("-f", "--orgfile", help="File of organism IDs to BLAST against (use this option if you want to test against multiple organisms", action="store", type="str", dest="orgfile", default=None)
parser.add_option("-d", "--db", help="BLAST database to use for BLASTing (use this option if you already have generated a blast database)", action="store", type="str", dest="db", default=None)
parser.add_option("-c", "--cutoff", help="E-value cutoff for TBLASTN (D=1E-5)", action="store", type="float", dest="cutoff", default=1E-5)
parser.add_option("-t", "--translation", help="Translation table number for TBLASTN (D=11 - bacteria, archaea and plant plastids)", action="store", type="int", dest="translation", default=11)
parser.add_option("-g", "--genecol", help="Column number for gene ID starting from 1 (D=1)", action="store", type="int", dest="gc", default=1)
parser.add_option("-r", "--orgcol", help="Column number for organism ID starting from 1 (D=1, ignored unless -f is specified)", action="store", type="int", dest="oc", default=1)
(options, args) = parser.parse_args()

#############
# Input sanity checks...
#############
if options.org is None and options.orgfile is None and options.database is None:
    sys.stderr.write("ERROR: One of -o (organism ID), -f (organism ID file), or -d (BLAST database) is a required input\n")
    exit(2)    

numnotnone = 0
if options.org is not None:
    numnotnone += 1
if options.orgfile is not None:
    numnotnone += 1
if options.db is not None:
    numnotnone += 1
if numnotnone > 1:
    sys.stderr.write("ERROR: Cannot specify more than one of -o, -f, and -d (input options)\n")
    exit(2)

#############
# Get a list of input protein IDs
#############
gc = options.gc - 1
pids = []
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    if spl[gc].startswith("fig|"):
        pids.append(spl[gc])
    else:
        pids.append("fig|%s" %(spl[gc]))

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Get protein sequences from the specified gene IDs
# and put them into a FASTA file for input into TBLASTN
rn = random.randint(0, 2**30)
qfile = "%d.faa" %(rn)
fid = open(qfile, "w")
q = "SELECT geneid, aaseq FROM processed WHERE geneid = ?;"
for pid in pids:
    cur.execute(q, (pid, ))
    for res in cur:
        fid.write(">%s\n%s\n" %(str(res[0]), str(res[1])))
fid.close()

# Compile a BLAST database to search against.
# The way we do this depends a bit on how the function was called.
# If an organism ID is passed... just find all the contigs from the database
# and compile them.
if options.org is not None:
    db = "%d.fna" %(rn)
    fid = open(db, "w")
    q = """SELECT contigs.contig_mod, contigs.seq FROM contigs
           WHERE contigs.organismid = ?"""
    cur.execute(q, (options.org, ) )
    atleastone = False
    for rec in cur:
        atleastone = True
        fid.write(">%s\n%s\n" %(rec[0], rec[1]))
    fid.close()
    if not atleastone:
        sys.stderr.write("ERROR: No contigs found for specified organism %s\n" %(options.org) )
        exit(2)
    sys.stderr.write("Making blast DB...\n")
    os.system("makeblastdb -dbtype nucl -in %s > /dev/null" %(db))
# If a file of organism IDs was specified, do the same thing but just do it for all of the specified IDs
elif options.orgfile is not None:
    db = "%d.fna" %(rn)
    oc = options.oc - 1
    fid_out = open(db, "w")
    q = """SELECT contigs.contig_mod, contigs.seq FROM contigs 
         WHERE contigs.organismid = ?"""
    for line in open(options.orgfile, "r"):
        spl = line.strip("\r\n").split("\t")
        orgid = spl[oc]
        cur.execute(q, (orgid, ) )
        atleastone = False
        for rec in cur:
            atleastone = True
            fid_out.write(">%s\n%s\n" %(rec[0], rec[1]))
        if not atleastone:
            sys.stderr.write("WARNING: No contigs found for specified organism %s\n" %(orgid))
            continue
    fid_out.close()
    # Compile BLAST database
    sys.stderr.write("Making blast DB...\n")
    os.system("makeblastdb -dbtype nucl -in %s > /dev/null" %(db))
elif options.db is not None:
    db = options.db

# Run TBLASTN
ofile = "%d.out" %(rn)
cmd = "tblastn -db %s -evalue %1.1e -query %s -out %s -db_gencode %d -outfmt \"6 qseqid sseqid sstart send evalue bitscore sframe\" " %(db, options.cutoff, qfile, ofile, options.translation)
sys.stderr.write("Now executing TBLASTN with command: \n%s\n" %(cmd))
os.system(cmd)

sys.stderr.write("Generating gene overlap report...\n")
q = "SELECT geneid, genestart, geneend, annotation, organism FROM processed WHERE contig_mod = ? AND MIN(genestart, geneend) <= MAX(?,?) AND MIN(?,?) <= MAX(genestart, geneend);"
for line in open(ofile, "r"):
    spl = line.strip("\r\n").split("\t")
    tblaststart = int(spl[2])
    tblastend = int(spl[3])
    subcontig = spl[1]
    cur.execute(q, (subcontig, tblaststart, tblastend, tblaststart, tblastend))
    atleastone = False
    for rec in cur:
        atleastone = True
        geneid = rec[0]
        genestart = int(rec[1])
        geneend = int(rec[2])
        annotation = rec[3]
        organism = str(rec[4])
        overlap = min(max(tblaststart, tblastend), max(genestart, geneend)) - max(min(tblaststart, tblastend), min(genestart, geneend))
        overlappct = float(overlap)/abs(float(genestart - geneend))*100
        # Same strand?
        if ( genestart < geneend and tblaststart < tblastend ) or ( genestart > geneend and tblaststart > tblastend ):
            print "%s\tSAMESTRAND\t%s\t%d\t%d\t%d\t%1.2f\t%s\t%s" %(line.strip("\r\n"), geneid, genestart, geneend, overlap, overlappct, annotation, organism)
        else:
            print "%s\tOTHERSTRAND\t%s\t%d\t%d\t%d\t%1.2f\t%s\t%s" %(line.strip("\r\n"), geneid, genestart, geneend, overlap, overlappct, annotation, organism)
    # No gene matches the loation
    if not atleastone:
        print "%s\tNOGENE\t\t\t\t\t\t\t" %(line.strip("\r\n"))

# Clean up
os.system("rm %d.*" %(rn))

con.close()
