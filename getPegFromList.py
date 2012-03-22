#!/usr/bin/python

# Get PEG from list of aliases since svr_alias_to_peg doesn't work but
# svr_aliases_of does.
import fileinput
import re

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    aliaslist = spl[1].split(",")
    for s in aliaslist:
        if re.match("fig\|\d+\.\d+", s):
            print "%s\t%s" % (spl[0], s)
