#!/usr/bin/python
#
# Use MAFFT to make a quick and dirty alignment between a (small) set of genes of interest.

import fileinput, optparse, os, random, sqlite3, sys
from locateDatabase import *

usage="%prog [options] < gene_list > alignment"
description="Get a quick and dirty (whole-gene) alignment between a set of genes of interest piped from stdin. Does not do any trimming."
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--genecol", help="Column number for gene starting from 1 (D=1)", action="store", type="int", dest="gc", default=1)
parser.add_option("-p", "--phylip", help="Output alignment as PHYLIP (D: FASTA)", action="store_true", dest="phylip", default=False)
parser.add_option("-k", "--keep", help="Keep temporary files (D: Delete them)", action="store_true", dest="keep", default=False)
(options, args) = parser.parse_args()

gc = options.gc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

id2seq = {}
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    gene = spl[gc]
    cur.execute("SELECT aaseq FROM processed WHERE processed.geneid = ?", (gene, ) )
    for res in cur:
        id2seq[gene] = str(res[0])

con.close()

# Make a temporary fasta file for MAFFT
rn = random.randint(0, 2**30)
fname = "%d.fasta" %(rn)
fid = open(fname, "w")
for myid in id2seq:
    fid.write(">%s\n%s\n" %(myid, id2seq[myid]))

fid.close()

# Run MAFFT
os.system("mafft --auto --reorder %s > %s.faa 2> /dev/null" %(fname, fname))

# If we want a PHYLIP file I'll give you one but it will have gene IDs replaced...
if options.phylip:
    sys.stderr.write("WARNING: Phylip file requested - due to need for 10-digit identifiers the gene IDs will be replaced.\n")
    os.system("cat %s.faa | fastaToPhylip.py -c %s.conv" %(fname, fname))
else:
    os.system("cat %s.faa" %(fname))

# Remove temporary files (including the alignment file)
if not options.keep:
    os.system("rm %d.*" %(rn) )
