#!/usr/bin/env python

# Generates a new table with only information we need to make it easier to calculate some statistics:
# Gene (col 1)
# Annotation (col 7)
# Also reads organism file (argument 1) and matches the strings.

import optparse
import sys

usage = "%prog Organism_file < RAW_file > Processed_file"
description="Processes organism file into a nicer format, adding organism and gene length information to the table and removing unneeded columns. The raw file must have organism / gene IDs similar to what is given by RAST and peg in teh appropriate column to designate proteins rather than nucleic acid entries."
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

if len(args) < 1:
    sys.stderr.write("ERROR: Organism file is a required input argument (see help text -h for details)\n")
    exit(2)

orgNames = {}
orgIds = {}

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
for line in open(args[0], "r"):
    spl = line.strip('\r\n').split('\t')
    # Organism ID is the key and orgnaism is the data
    # Critical that we include the extra "." here because we want to 
    # search for only the ID itself, which is always followed by a . in the SEED IDs
    # For example, fig|333.3.peg.3
    # 333.3 is the ID and we don't want to match 333.33 or 3333.3
    orgNames["fig|" + spl[1] + "."] = spl[0]
    orgIds["fig|" + spl[1] + "."] = spl[1]

import fileinput
for line in fileinput.input("-"):
    spl = line.strip('\r\t').split('\t')
    # Protein-coding DNAs only (this also cuts off the title row)
    if not spl[2] == "peg":
        continue

    # Search for which organism [this is probably always the same]
    myorg = ""
    myid = ""
    
    for org in list(orgNames.keys()):
        if org in spl[1]:
            myorg = orgNames[org]
            myid = orgIds[org]
            break

    contig = myid + "." + spl[0]

    modln = spl[1] + "\t" + myid + "\t" + myorg + "\t\t" + contig + "\t" + str(sign(int(spl[5]) - int(spl[4]))) + "\t" + str(len(spl[11])) + "\t" + str(len(spl[12]))
    print(modln)

