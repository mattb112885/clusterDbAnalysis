#!/usr/bin/python
#
# This is  a pipe command
# Create a presence-absence fasta file
# from teh presence-absence table
# The presence-absence table starts with interesting stuff on the third column

import fileinput

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

print orgdict

for key in presabsdict:
    print ">%s" %(orgdict[key])
    print "".join(str(s) for s in presabsdict[key])
