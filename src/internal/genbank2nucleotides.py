#!/usr/bin/python

import fileinput, optparse, re, sys
from Bio import SeqIO

usage="%prog -f [genbank_file] > fna file"
description="""
Make a contig nucleic acid FASTA file out of a genbank file. Contig names MUST have no spaces.
You should pass to this an organism ID so that the contig names are unique for particular organisms (and so that they match what is in the database)
If no organism ID is passed and it cannot be inferred from the file name we throw an error.
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-t", "--tab", help="Instead of a FASTA file, print a tab-delimited file with contig in column 1 and sequence in column 2", action="store_true", dest="tab", default=False)
parser.add_option("-o", "--org", help="Organism ID (e.g. 83333.1) (D=Try to read from filename)", action="store", type="str", dest="orgid", default=None)
parser.add_option("-f", "--file", help="genbank file [D: None]", action="store", type="str", dest="genbank", default=None)
(opts, args) = parser.parse_args()

if opts.genbank == None:
    sys.stderr.write("ERROR: Genbank file (-f) must be provided\n")
    exit(2)

orgname = ""
if opts.orgid is None:
    sys.stderr.write("WARNING: No organism ID provided. Attempting to read the organism ID from the filename\n")
    idFinder = re.compile("\d+\.\d+")
    mtch = idFinder.search(opts.genbank)
    if mtch is None:
        # Organism name must be specified in order to import any of this info into the database. Since that
        # is the primary purpose of this file, we make this an error rather than a warning
        sys.stderr.write("ERROR: No IDs in expected format \d+\.\d+ found in file name. You must specify -o to append organism IDs to the contig names\n")
        exit(2)
    else:
        orgname = mtch.group(0)
        sys.stderr.write("Found organism ID: %s\n" %(orgname))
else:
    orgname = opts.orgid

# Note - biopython can handle multi-genbanks...
multi_records = SeqIO.parse(opts.genbank, "genbank")

for gb_seqrec in multi_records:
    # Note - this MUST be consistent with convertGenbankToTable.py or else the contig names won't match the ones from the raw files and you won't be able to do anything with the loaded sequences.
    if gb_seqrec.id == "unknown":
        contig_name = gb_seqrec.name
    else:
        contig_name = gb_seqrec.id
    contig_name = "%s.%s" %(orgname, contig_name)
    seq = str(gb_seqrec.seq)
    if opts.tab:
        print "%s\t%s\t%s" %(contig_name, seq, orgname)
    else:
        print ">%s\n%s" %(contig_name, gb_seqrec.seq)
