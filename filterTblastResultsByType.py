#!/usr/bin/python

# This is a pipe command - pipe in the results from the db_clusterToTblastn.py
#
# Line by line, sorts the results by their type
# Type is determined by two ratios:
#   R1: TBLASTN length / Query length
#   R2: Target length / TBLASTN length
#
# TYPES:
# SAME - lb <= R1 <= ub
#        lb <= R2 <= ub
# TFRAGQ - R1 <= lb
#         lb <= R2 <= ub
# DOMAIN - R1 <= lb
#          ub <= R2
# QFRAGT - lb <= R1 <= ub
#          ub <= R2
# EARLYSTOP - lb <= R1 <= ub
#             R2 <= lb
# COMBINATION: R1 <= lb
#              R2 <= lb
#
# SAME: Predict that it is not a pseudogene (the two genes are the same)
# TFRAGQ: Predict that target is a fragment of the query
# DOMAIN: Predict that the TBLASTN hit was only a domain-specific hit.
# QFRAGT: Predict that the query is a fragment of the target
# EARLYSTOP: Predict an early STOP codon in the target gene
# COMBINATION: Some uncharacterized combination
#
# Does not print SAME - only prints the others to save space
# and make it easier to aggregate statistics

lb = 0.75
ub = 1.25

import fileinput

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    # abs is needed because negatives occur due to strand differences...
    R1 = abs(float(spl[13])/float(spl[5]))
    R2 = abs(float(spl[16])/float(spl[13]))

    typestring = ""

    if lb <= R1 and R1 <= ub and lb <= R2 and R2 <= ub:
        # Do not print SAME strings?
        typestring = "SAME"
        # continue
    elif R1 <= lb and lb <= R2 and R2 <= ub:
        typestring = "TFRAGQ"
    elif R1 <= lb and ub <= R2:
        typestring = "DOMAIN"
    elif lb <= R1 and R1 <= ub and ub <= R2:
        typestring = "QFRAGT"
    elif lb <= R1 and R1 <= ub and R2 <= lb:
        typestring = "EARLYSTOP"
    elif R1 <= lb and R2 <= lb:
        typstring = "COMBINATION"
    else:
        typestring = "??"

    print "\t".join(spl) + "\t" + str(R1) + "\t" + str(R2) + "\t" + typestring
