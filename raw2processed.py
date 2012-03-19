#!/usr/bin/python

# Generates a new table with only information we need to make it easier to calculate some statistics:
# Gene (col 1)
# Annotation (col 7)
# Also reads organism file (argument 1) and matches the strings.

import sys

fname = sys.argv[1]
fid = open(fname, "r")

orgNames = {}
orgIds = {}
orgAbbrevs = {}

def sign(num):
    if num == 0:
        return 0
    elif num < 0:
        return -1
    elif num > 0:
        return 1
    else: #NAN / None
        return None

# Organisms file
for line in fid:
    spl = line.strip().split('\t')
    # Organism ID is the key and orgnaism is the data
    # Critical that we include the extra "." here because we want to 
    # search for only the ID itself, which is always followed by a . in the SEED IDs
    # For example, fig|333.3.peg.3
    # 333.3 is the ID and we don't want to match 333.33 or 3333.3
    orgNames["fig|" + spl[2] + "."] = spl[0]
    orgIds["fig|" + spl[2] + "."] = spl[2]
    orgAbbrevs["fig|" + spl[2] + "."] = spl[1]

import fileinput
for line in fileinput.input("-"):
    spl = line.strip().split('\t')
    # Protein-coding DNAs only (this also cuts off the title row)
    if not spl[2] == "peg":
        continue

    # Search for which organism [this is probably always the same]
    myorg = ""
    myid = ""
    myabbrev = ""
    
    for org in orgNames.keys():
        if org in spl[1]:
            myorg = orgNames[org]
            myid = orgIds[org]
            myabbrev = orgAbbrevs[org]
            break

    contig = myid + "." + spl[0]

    modln = spl[1] + "\t" + myid + "\t" + myorg + "\t" + myabbrev + "\t" + contig + "\t" + str(sign(int(spl[5]) - int(spl[4]))) + "\t" + str(len(spl[11])) + "\t" + str(len(spl[12]))
    print modln

fid.close()
