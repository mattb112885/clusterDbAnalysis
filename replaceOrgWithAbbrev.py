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

if len(sys.argv) < 2:
    print "Usage: ./replaceOrgWithName.py [orgfile] [useAbbrev]"
    print "./replaceOrgWithName [orgfile] [TRUE] means use the organism abbreviation rather than the full name"
    print "Omit useAbbrev to make it FALSE (or use anything but TRUE or T)"
    exit(2)

useabbrev = False
if len(sys.argv) == 3:
    if sys.argv[2].lower().startswith("t"):
        useabbrev = True

orgAbbrev = {}
fid = open(sys.argv[1], "r")
for line in fid:
    spl = line.strip().split("\t")
    if useabbrev:
        orgAbbrev[spl[2]] = spl[1].replace(" ", "_")
    else:
        orgAbbrev[spl[2]] = spl[0].replace(" ", "_")

for line in fileinput.input("-"):
    myline = line.strip().replace("fig|", "")
    for s in orgAbbrev:
        myline = myline.replace(s, orgAbbrev[s])
    print myline
