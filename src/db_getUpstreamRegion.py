#!/usr/bin/env python

# From my initial tests I think the results are correct but it's still possible
# an off-by-one is lurking somewhere...

import fileinput, optparse, sqlite3, sys
from FileLocator import *
# For sequence transposing
from Bio import Seq

header = [ "geneid", "status", "upstream_sequence" ]
usage = """%prog [options] < geneids > geneids_with_upstream

Output table: """ + " ".join(header)
description="""Get the upstream nucleotide sequence of the given set of genes.
Requires you to have the contigs loaded into the database. It is careful to only
print sequences up to the next called gene or gap (N) in the sequence unless you
tell it not to. This is because upstream regions are often used to search for motifs
but motifs may be hard to distinguish from strong protein sequence conservation in
protein-coding regions.

Explanation of warnings:
CONTIGEND - Upstream region reached the end of a contig.
NOUPSTREAM - There is no upstream region (either another gene overlaps or immediately abuts, or the gene reaches the very end of a contig)
CONTAINSGAP - if there are a large number of n's (user-defined but default is to warn if any N's are present).
OTHERGENE - If there is another gene within N nucleotides upstream of the query gene.

Automatically transposes the DNA and rotates it around so that it should be possible to directly compare motifs...
"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--numupstream", help="Number of NT to go upstream (D = 100). For 'up to the next gene' use a large number and do not specify -o", action="store", type="int", dest="numupstream", default=100)
parser.add_option("-w", "--gapwarn", help="Number of N's in the upstream region before warning about a gap (D = 1, warn if any Ns are present)", action="store", type="int", dest="gapwarn", default=1)
parser.add_option("-o", "--allowgeneoverlap", 
                  help="If specified, always try to reach the number of upstream nucleotides even if another called gene is there (D: Cut it off with warning)", 
                  action="store_true", dest="allowgeneoverlap", default=False)
parser.add_option("-l", "--othergenelength", help="If allowing gene overlaps, still ignore called genes less than this length (in nucleotides) within the upstream region (D=0 - cut off after ANY gene)", action="store", type="int", dest="othergenelength", default=0)
parser.add_option("-g", "--genecolumn", help="Column number for gene ID in input, starting from 1 (D=1)", action="store", type="int", dest="gc", default=1)
parser.add_option("-i", "--ingene", help="Number of nucleotides WITHIN the gene to gather in addition to the upstream region (D=3 - i.e. grab the start codon only)", action="store", type="int", dest="ingene", default=3)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

# Read stdin to get gene IDs
gc = options.gc - 1
geneids = []
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    if spl[gc].startswith("fig|"):
        geneids.append(spl[gc])
    else:
        geneids.append("fig|%s" %(spl[gc]) )

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Lets get a dict from contig ID to contig - this saves us having to do JOINs many times with long strings.
contigToSeq = {}
q = "SELECT contig_mod, seq FROM contigs;"
cur.execute(q)
for rec in cur:
    contigToSeq[rec[0]] = rec[1]

q1 = "SELECT geneid, genestart, geneend, strand, processed.contig_mod FROM processed WHERE geneid = ?"
# The ? here are the start and stop sites of the upstream region (which must not overlap with the query gene)
#, the contig for the query gene, and the specified minimum nucleotide length for the called ORF
q2 = """SELECT geneid, genestart, geneend, strand FROM processed 
        WHERE MAX(genestart, geneend) >= MIN(?,?) AND MIN(genestart, geneend) <= MAX(?,?) AND contig_mod = ?
        AND ABS(genestart - geneend) > ?"""
NUMINGENE = options.ingene

if options.header:
    print "\t".join(header)
for gene in geneids:
    # Where is this gene located?
    cur.execute(q1, (gene, ) )
    warnings = ""
    # I use fetchall() here because I need to re-use the cursor inside the loop.
    for rec in cur.fetchall():
        if not rec[4] in contigToSeq:
            sys.stderr.write("WARNING: The contig for gene %s was not found in the database so it will be skipped\n" %(gene))
            continue
        contigseq = contigToSeq[rec[4]]
        contiglen = len(contigseq)
        # Get the start and stop locations for the desired upstream region
        strand = rec[3]
        if strand == "-":
            # Start > stop for - strand genes
            interval_start = int(rec[1]) + 1
            interval_end = interval_start + options.numupstream - 1
        else:
            interval_end = int(rec[1]) - 1
            interval_start = interval_end - options.numupstream + 1
        # This means the gene is at the very end of the contig so lets just kill it.
        if interval_end < 1 or interval_start > contiglen:
            warnings += "NOUPSTREAM,"
            print "%s\t%s\t" %(geneid,warnings)
            continue
        if interval_start < 1 or interval_end > contiglen:
            warnings += "CONTIGEND,"
            if interval_start < 1:
                interval_start = 1
            if interval_end > contiglen:
                interval_end = contiglen

        # Is there a called gene in the interval? (with at least the specified minimum length)
        cur.execute(q2, (interval_start, interval_end, interval_start, interval_end, rec[4], options.othergenelength))
        if strand == "-":
            interval_start = interval_start - NUMINGENE
        else:
            interval_end = interval_end + NUMINGENE

        s = cur.fetchall()
        if len(s) > 0:
            warnings += "OTHERGENE,"
            if not options.allowgeneoverlap:
                for r in s:
                    if strand == "-":
                        if min(int(r[1]), int(r[2])) < interval_start:
                            warnings += "NOUPSTREAM,"
                            interval_end = interval_start - 1
                        else:
                            # In case there are multiple called genes in the interval, and since we replaced
                            # interval_end in a previous iteration, we need to make sure we allow
                            # the replaced one to still be the lowest.
                            interval_end = min(interval_end, int(r[1])-1)
                            interval_end = min(interval_end, int(r[2])-1)
                    else:
                        if max(int(r[1]), int(r[2])) > interval_end:
                            warnings += "NOUPSTREAM,"
                            interval_start = interval_end + 1
                        else:
                            interval_start = max(interval_start, int(r[1])+1)
                            interval_start = max(interval_start, int(r[2])+1)
        # Get the actual sequence from what is left.
        # BUT: Note that the start and stop locations start at 1, while the
        # array indexes start at 0!
        startidx = interval_start - 1
        stopidx = interval_end - 1
        seq = contigseq[startidx:stopidx+1]
        if seq.lower().count("n") > options.gapwarn:
            warnings += "CONTAINSGAP"
        # We need to do the reverse complement if we are on the "-" strand...
        if strand == "-":
            seq = str(Seq.Seq(seq).reverse_complement())
        print "%s\t%s\t%s" %(gene, warnings, seq)

con.close()
