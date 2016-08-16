#!/usr/bin/env python

from __future__ import print_function
import optparse
import sqlite3
import sys

from FileLocator import *
from ClusterFuncs import *

usage = "%prog [options] > contig_ids"
description = """Get contig IDs. 
By default returns ALL contig IDs. Optionally return contigs only
for specific organisms.
"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-i", "--organismid", help="Only return contigs for organisms matching this organism ID", action="store", type="str", dest="organismid", default=None)
parser.add_option("-o", "--organism", help="Only return contigs with the specified organism name", action="store", type="str", dest="organism", default=None)
parser.add_option("-s", "--sanitized", help="Specify this if the organism's name is sanitized", action="store", type="str", dest="sanitized", default=False)
(options, args) = parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

if options.organismid is not None and options.sanitized:
    raise NotImplementedError("I haven't yet implemented putting a sanitized organism ID as an input")

contiglist = getContigIds(cur, orgid=options.organismid, orgname=options.organism, issanitized=options.sanitized)

for contig in contiglist:
    print(contig)

con.close()
