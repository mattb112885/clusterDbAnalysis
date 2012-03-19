#!/usr/bin/python

# Default output is to stdout
# NOTE - this function will fail if any of the contig names are too long because Biopython
# complains that the locus tags are malformed and then it crashes.
# 101192.3 works OK

from Bio import SeqIO
import optparse
import sys

usage = "%prog -gb [genbankefile] -tsv [tab-delimited file from RAST]"
description = "Given list of genes to match, returns a list of BLAST results between genes in the list only"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("--tsv", help="Tab-delimited file from RAST", action="store", type="string", dest="tsv")
parser.add_option("--gb",  help="GENBANK file from RAST", action="store", type="string", dest="genbank")
parser.add_option("--o", help="Output file (default = stdout)", action="store", type="string", dest="outfile", default=None)

(options, args) = parser.parse_args()

if options.tsv == None or options.genbank == None:
    print "ERROR: Must provide both SEED tsv file and genbank file"
    exit(2)

tsv = open(options.tsv, "r")
# Gene ID is on the 2nd and AA on the 13th column.
seqs = {}
for line in tsv:
    spl = line.strip().split("\t")
    # Skip over things with no AA sequence like RNA genes...
    if(len(spl) < 13):
        continue
    # Use sequence as the key
    seqs[spl[12]] = spl[1]
tsv.close()

# Add gene IDs to the parsed genbank file based on the AA sequences...
gbk = open(options.genbank, "r")
records = list(SeqIO.parse(gbk, "genbank"))

for s in records:
    for l in s.features:
        if l.type == "CDS":
            seq = l.qualifiers["translation"][0]
            # Print to make sure i'm donig this right...
            if seq in seqs:
                l.qualifiers["locus_tag"] = [seqs[seq]]

if options.outfile == None:
    outhandle = sys.stdout
else:
    outhandle = open(options.outfile, "w")

SeqIO.write(records, outhandle, "genbank")
