#!/usr/bin/python

# Get PEG from list of aliases since svr_alias_to_peg doesn't work but
# svr_aliases_of does.
import fileinput
import re

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    # If no match is found we don't just want to get an error message...
    if len(spl) < 3:
#        print "%s\t%s\t%s" %(spl[0], spl[1], "")
        continue

    aliaslist = spl[2].split(",")
    for s in aliaslist:
        if re.match("fig\|\d+\.\d+", s):
            print "%s\t%s\t%s" % (spl[0], spl[1], s)
