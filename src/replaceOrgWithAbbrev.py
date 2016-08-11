#!/usr/bin/env python

# Usage: replaceOrgWithAbbrev.py [orgfile]
#
# Pipe in a text file and it will replace the organism ID (fig|xx.yy.zz) with
# the organism abbreviation in all cases where it appears.
#
# Results are printed to stdout.
# Example usage:
# - Pipe in a Newick file containing gene IDs (or just organism IDs)
# and it will make them more readable by replacing them with organism names
#

import fileinput
import sys
import optparse
from sanitizeString import *
from FileLocator import *

usage="%prog [options] < text_file > text_file_with_orgname"
description="Replace organism IDs (fig|xx.yy) with organism abbreviations in a text file (e.g. a newick file). Note the fig| is optional."
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--orgfile", help="Organism file (optional, default = organisms file in root directory of current install)", 
                  action="store", type="str", dest="orgfile", default=None)
parser.add_option("-k", "--keeppeg", help="Keep PEG ID? (if specified, keeps peg id. If not, throws it away)", action="store_true", dest="keeppeg", default=False)
parser.add_option("-s", "--sanitized", help="Specify this if the organism IDs are sanitized in the file (fig_xx_yy instead of fig|xx.yy)", 
                  action="store_true", dest="sanitized", default=False)
(options, args) = parser.parse_args()

if options.orgfile is None:
    options.orgfile = locateOrganismFile()

keeppeg = options.keeppeg

orgAbbrev = {}
fid = open(options.orgfile, "r")
for line in fid:
    spl = line.strip('\r\n').split("\t")
    if options.sanitized:
        orgid = sanitizeString(spl[1], False)
    else:
        orgid = spl[1]

    orgAbbrev[orgid] = sanitizeString(spl[0], False)

for line in fileinput.input("-"):
    myline = line.strip('\r\n')
    # I always replace this since it shouldn't break anything to leave it out anyway
    myline = myline.replace("fig|", "")

    for s in orgAbbrev:
        if keeppeg:
            myline = myline.replace(s, orgAbbrev[s] + "_" + s)
        else:
            myline = myline.replace(s, orgAbbrev[s])
    print(myline)
