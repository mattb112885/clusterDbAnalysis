#!/usr/bin/python

# Wrapper script for TBLASTN verification of presence\absence
# across different genomes

import fileinput, os, optparse, random, sqlite3, sys
from locateDatabase import *

usage = "%prog -d Fasta_db [options] < Protein_ids > Tblastn_table"
description = """ Given a Fasta DB (contig names in that fasta file must be
the SAME as the contig names in the database), and a set of gene IDs, 
attempts to run TBLASTN and identify
missing genes. It identifies called genes that match the hit location and
also tries to find genes on the opposite strand that conflict."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-d", "--db", help="BLAST database for running TBLASTN (mandatory)", action="store", type="str", dest="db", default=None)
parser.add_option("-c", "--cutoff", help="E-value cutoff for TBLASTN (D=1E-5)", action="store", type="float", dest="cutoff", default=1E-5)
parser.add_option("-t", "--translation", help="Translation table number for TBLASTN (D=11 - bacteria, archaea and plant plastids)", action="store", type="int", dest="translation", default=11)
parser.add_option("-g", "--genecol", help="Column number for gene ID starting from 1 (D=1)", action="store", type="int", dest="gc", default=1)
(options, args) = parser.parse_args()

#############
# Set up database call for input proteins
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
rn = random.randint(0, 2**30)
qfile = "%d.faa" %(rn)
fid = open(qfile, "w")
q = "SELECT geneid, aaseq FROM processed WHERE geneid = ?;"
for pid in pids:
    cur.execute(q, (pid, ))
    for res in cur:
        fid.write(">%s\n%s\n" %(str(res[0]), str(res[1])))
fid.close()

# Run TBLASTN
ofile = "%d.out" %(rn)
cmd = "tblastn -db %s -evalue %1.1e -query %s -out %s -db_gencode %d -outfmt \"6 qseqid sseqid sstart send evalue bitscore sframe\" " %(options.db, options.cutoff, qfile, ofile, options.translation)
sys.stderr.write("Now executing TBLASTN with command: \n%s\n" %(cmd))
os.system(cmd)

sys.stderr.write("Generating gene overlap report...\n")
q = "SELECT geneid, genestart, geneend FROM processed WHERE contig_mod = ? AND MIN(genestart, geneend) <= MAX(?,?) AND MIN(?,?) <= MAX(genestart, geneend);"
for line in open(ofile, "r"):
    spl = line.strip("\r\n").split("\t")
    tblaststart = int(spl[2])
    tblastend = int(spl[3])
    subcontig = spl[1]
    print subcontig, tblaststart, tblastend
    cur.execute(q, (subcontig, tblaststart, tblastend, tblaststart, tblastend))
    atleastone = False
    for rec in cur:
        atleastone = True
        geneid = rec[0]
        genestart = int(rec[1])
        geneend = int(rec[2])
        overlap = min(max(tblaststart, tblastend), max(genestart, geneend)) - max(min(tblaststart, tblastend), min(genestart, geneend))
        overlappct = float(overlap)/abs(float(genestart - geneend))
        # Same strand?
        if ( genestart < geneend and tblaststart < tblastend ) or ( genestart > geneend and tblaststart > tblastend ):
            print "%s\tSAMESTRAND\t%s\t%d\t%d\t%d\t%1.2f" %(line.strip("\r\n"), geneid, genestart, geneend, overlap, overlappct)
        else:
            print "%s\tOTHERSTRAND\t%s\t%d\t%d\t%d\t%1.2f" %(line.strip("\r\n"), geneid, genestart, geneend, overlap, overlappct)
    # No gene matches the loation
    if not atleastone:
        print "%s\tNOGENE\t\t\t\t\t" %(line.strip("\r\n"))

# Clean up
os.system("rm %d.*" %(rn))

con.close()
