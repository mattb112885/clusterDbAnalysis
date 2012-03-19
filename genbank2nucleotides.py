#!/usr/bin/python

# This is a pipe command. Pipe in a genbank file and
# get out a nucleotide sequence FASTA file containing
# all the contigs in the genome

# Requires biopython.
# Note - the actual stuff we want is printed to stdout.
#   Python will generate some errors because the seed exports "UNK" as the date of
#   modification (in the LOCUS field) instead of an actual valid date.
#   As far as I can tell, these can be safely ignored.

from Bio import SeqIO
import sys

records = SeqIO.parse(sys.stdin, "genbank")

# Make a nucleotide FASTA file and export to stdout.
# Make sure the genome ID is exported... with the extra
# period to make sure we don't only match subsets of the ID
genome_id = ""
for s in records:
    if genome_id == "":
        for l in s.features:
            if l.type == "source":
                genome_id = l.qualifiers["genome_id"][0]
                break
    print ">" + genome_id + "." + s.name
    print s.seq
