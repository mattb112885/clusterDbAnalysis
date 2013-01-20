#!/usr/bin/python

import fileinput, optparse, sqlite3, sys
from FileLocator import *
from getSequenceRegion import *

usage = "%prog [options] < table > table_with_sequences"
description = """Identify the sequence of the HSP based on
the provided ID (contig or gene) and the start and stop locations
of the HSP.

The user must specify one of -p, -n, and -c - it gets too ambiguous
otherwise (especially between -p and -n). This tells the program what type
of data to use to obtain the result.

-p is for BLASTP results, -n is for BLASTN results against genes,
and -c is for getting the sequence for a certain location on a contig from
BLASTN\TBLASTN results.

For DNA sequences, start < stop is assumed to mean + strand and
stop < start means - strand.

If the "translate DNA" (-t) option is specified, it will pull out the DNA
sequence, do the reverse-complement if necessary, and then translate that.
"""

parser = optparse.OptionParser(usage=usage, description=description)

# Input options
parser.add_option("-n", "--nucleotide", help="Specify this flag if the piped-in locations are on _genes_ (i.e. BLASTN vs genes rather than BLASTP against proteins)", action="store_true", dest="nucl", default=False)
parser.add_option("-c", "--contig", help="Specify this flag if the piped-in locations are on a _contig_ (i.e. from BLASTN or TBLASTN vs the whole contig)", action="store_true", dest="contig", default=False)
parser.add_option("-p", "--protein", help="Specify this flag is the piped-in locations are within a _protein_ sequence", action="store_true", dest="prot", default=False)

parser.add_option("-i", "--idcol", help="Column number for ID to use starting from 1 (D: 2 as from BLASTP or BLASTN - use 3 for TBLASTN)", action="store", type="int", dest="idcol", default=2)
parser.add_option("-s", "--startcol", help="Column number for start location to use - NOTE this will be different from the default for TBLASTN results (D: 9 - target start location for BLASTP and BLASTN. Use 5 for TBLASTN)", action="store", type="int", dest="startcol", default=9)
parser.add_option("-e", "--endcol", help="Column number for end \ stop location to use - NOTE this will be different from the default for TBLASTN results (D: 10 - target stop location for BLASTP and BLASTN. Use 6 for TBLASTN)", action="store", type="int", dest="endcol", default=10)

# Output options
parser.add_option("-t", "--translate_dna", help="Given a set of options that would result in a DNA sequence, translate that sequence and return the translated sequence instead. (D: If relevant, just returns the DNA sequence) - WARNING: Do not use for BLASTN since it can cause issues with translation frame! Use for TBLASTN only", action="store_true", dest="translate", default=False)

(options, args) = parser.parse_args()

if options.prot and options.translate:
    sys.stderr.write("WARNING: With protein locations as input the -t option does not make sense - it will be ignored\n")
    options.translate = False

############################
# Get requested sequences  #
############################

# Construct the query to get the starting sequence
if options.contig:
    q = "SELECT seq FROM contigs WHERE contig_mod = ?"
    isprot = False
elif options.nucl:
    q = "SELECT nucseq FROM processed WHERE geneid = ?"
    isprot = False
elif options.prot:
    q = "SELECT aaseq FROM processed WHERE geneid = ?"
    isprot = True
else:
    sys.stderr.write("ERROR: One of -c, -p, and -n must be specified. Type -h for details\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

idcol = options.idcol - 1
startcol = options.startcol - 1
endcol = options.endcol - 1

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    try:
        myid = spl[idcol]
        sloc = int(spl[startcol])
        eloc = int(spl[endcol])
    except IndexError:
        sys.stderr.write("ERROR: Input table does not have enough columns for the specified column numbers for ID, start and stop location\n")
        exit(2)
    except ValueError:
        sys.stderr.write("ERROR: The specified column numbers for start and stop locations do not appear to have integers in them - check that the specified column numbers are correct\n")
        exit(2)

    cur.execute(q, (myid,))
    for rec in cur:
        startingSeq = str(rec[0])
        finalSeq = getContigRegion(startingSeq, sloc, eloc, isprot, options.translate)
        print "\t".join(spl + [ finalSeq ] )
