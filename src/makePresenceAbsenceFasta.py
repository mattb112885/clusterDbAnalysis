#!/usr/bin/env python
#
# This is  a pipe command
# Create a presence-absence fasta file
# from teh presence-absence table
# The presence-absence table starts with interesting stuff on the third column

import fileinput
import optparse

usage = "%prog < presence_absence_01_table > presence_absence_fasta_file"
description="Create a presence/absence fasta file from a presence / absence 0-1 table"
parser = optparse.OptionParser(usage=usage, description=description)
parser.parse_args()

# On the first line lie the labels.
firstline = True

presabsdict = {}
orgdict = {}

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    if firstline:
        for s in range(3,len(spl)):
            presabsdict[s] = []
            orgdict[s] = spl[s]
        firstline = False
        continue

    for s in range(3, len(spl)):
        presabsdict[s].append(int(spl[s]))

#print orgdict

for key in presabsdict:
    print(">%s" %(orgdict[key]))
    print("".join(str(s) for s in presabsdict[key]))
