#!/usr/bin/python

import fileinput, optparse, re, sys
from sanitizeString import *

usage="%prog [ID] < [Infile] > Raw_file"
description="""This is a function for internal use. Call convertKbaseToRast.sh instead.
This function takes an intermediate file created by convertKbaseToRast.sh and
an assigned genome ID, and converts all of the IDs to the approrpiate format for use
with the database.
"""

parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

organismid = args[0]
orgpattern = re.compile("\d+\.\d+")
if orgpattern.search(organismid) is None:
    sys.stderr.write("ERROR: Specified organism ID %s does not match the expected format for organism IDs\n" %(organismid))
    exit(2)

pattern = re.compile("(peg|CDS)\.(\d+)")

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    ftype = spl[7]
    # We only suppoer protein-coding regions.
    # CDS and peg are synonyms in the KBase.
    if not (ftype == "CDS" or ftype == "peg"):
        continue
    ftype = "peg"
    contig_id = spl[1].replace("kb|", "")
    annotation = spl[6] 
    geneid = spl[5]
    # Add gene ID to the annotation for easier back-searching
    annotation = annotation + "_" + geneid
    # Convert gene ID to have format fig|organismid.peg.\d+
    match = pattern.search(geneid)
    # Just in case.
    if match is None:
        continue
    geneid = "fig|%s.peg.%s" %(organismid, match.group(2))

    # Deal with the start and stop locations
    # The KBase standard is different from the RAST one - the "start"
    # is really the first base on the contig.
    #
    firstbase = int(spl[2])
    length = int(spl[3])
    strand = spl[4]
    if strand == "+":
        start = firstbase
        stop = firstbase + length - 1
    elif strand == "-":
        stop = start - length + 1

    # Get other stuff needed for the RAW file
    dnaseq = spl[8]
    aaseq = spl[9]

    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(contig_id, geneid, "peg", contig_id, start, stop, strand, annotation, "", "", "", dnaseq, aaseq)
