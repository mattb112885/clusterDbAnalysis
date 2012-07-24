#!/usr/bin/python

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

usage="%prog -f orgfile [options] < text_file > text_file_with_orgname"
description="Replace organism IDs (fig|xx.yy) with organism abbreviations in a text file (e.g. a newick file)"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--orgfile", help="Organism file (required)", action="store", type="str", dest="orgfile", default=None)
parser.add_option("-a", "--useabbrev", help="Use abbreviation? (If specified, use the abbreviated form of the organism name. If not specified, use the entire organism name)", 
                  action="store_true", dest="useabbrev", default=False)
parser.add_option("-k", "--keeppeg", help="Keep PEG ID? (if specified, keeps peg id. If not, throws it away)", action="store_true", dest="keeppeg", default=False)
(options, args) = parser.parse_args()

if options.orgfile == None:
    sys.stderr.write("ERROR: Orgfile (-f orgfile) is a required argument to replaceOrgWithAbbrev\n")
    exit(2)

useabbrev = options.useabbrev
keeppeg = options.keeppeg

orgAbbrev = {}
fid = open(options.orgfile, "r")
for line in fid:
    spl = line.strip('\r\n').split("\t")
    if useabbrev:
        orgAbbrev[spl[2]] = spl[1].replace(" ", "_")
    else:
        orgAbbrev[spl[2]] = spl[0].replace(" ", "_")

for line in fileinput.input("-"):
    myline = line.strip('\r\n')

    if not keeppeg:
        myline = myline.replace("fig|", "")

    for s in orgAbbrev:
        if keeppeg:
            myline = myline.replace(s, orgAbbrev[s] + "_" + s)
        else:
            myline = myline.replace(s, orgAbbrev[s])
    print myline
