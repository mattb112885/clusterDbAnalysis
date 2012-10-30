#!/usr/bin/python

import optparse
import os
import re
import sys

usage="%prog [organism_file_with_abbrevs] [organism_file_to_replace_abbrevs]"
description = """This function is intended for internal use.
Merge abbreviations from the old organism file into the new one
in cases where they exist in the old one. 
Replaces the new organism file (without warning) with another one
hat has those abbreviations in place."""
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

if len(args) != 2:
    sys.stderr.write("ERROR: both old and new organism files are required\n")
    exit(2)

oldorgfile = args[0]
neworgfile = args[1]

if oldorgfile == neworgfile:
    sys.stderr.write("ERROR: Old organism file must be different from new organism file.\n")
    exit(2)

org2abbrev = {}
for line in open(oldorgfile, "r"):
    spl = line.strip("\r\n").split("\t")
    # We do NOT want to merge default identifiers (DEFAULT_##)
    # because then we can get conflicts as new organisms
    # are added
    if re.match("^DEFAULT_\d+$", spl[1]) is not None:
        sys.stderr.write("WARNING: Ignoring default abbreviation %s in the old organism file\n" %(spl[1]))
        continue
    org2abbrev[spl[0]] = spl[1]

tmporgfile = "TEMPORARY_ORGANISM_FILE"
tmporgfid = open(tmporgfile, "w")

for line in open(neworgfile, "r"):
    spl = line.strip("\r\n").split("\t")
    if spl[0] in org2abbrev:
        spl[1] = org2abbrev[spl[0]]
    tmporgfid.write("\t".join(spl) + "\n")

tmporgfid.close()
os.system("mv %s %s" %(tmporgfile, neworgfile))
