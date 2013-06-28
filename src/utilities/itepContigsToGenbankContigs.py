#!/usr/bin/python

import fileinput
import optparse
import re
import sys

usage = "%prog (options) < [ITEP_contig_ids] > Genbank_contig_ids"
description = """Convert a list of ITEP contig IDs to a list of contig IDs
that would match them in the Genbank file for the same organism.
We recommend ONLY running this on one organism's contigs at a time due to
possibility of namespace collision (multiple organisms could have teh same
contig IDs in the genbank files)."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-t", "--truncate", help="Specify this flag if the Genbank file has truncated (to 16 characters) contig IDs BUT ITEP does not",
                  action="store_true", dest="truncate", default=False)
parser.add_option("-c", "--contig_col", help="Column number for contig starting from 1 (D:1)", action="store", type="int", dest="cc", default=1)
(options, args) = parser.parse_args()

cc = options.cc - 1

organismIdPlusPeriod = re.compile("^\d+\.\d+\.")

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    itep_contig = spl[cc]
    # ITEP contig IDs are the Genbank contig IDs appended to the organism ID.
    # organismID.contigID
    # Strip off the organism ID
    raw_contig = organismIdPlusPeriod.sub("", itep_contig)
    if options.truncate:
        genbank_contig = raw_contig[:16]
    else:
        genbank_contig = raw_contig
        pass
    print "%s\t%s" %(itep_contig, genbank_contig)
