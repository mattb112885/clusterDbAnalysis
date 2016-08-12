#!/usr/bin/env python

from __future__ import print_function
import re

def convertGprListToGeneList(rxn2gpr, searchkey="fig\|\d+\.\d+\.peg\.\d+"):
    '''Convert a set of gene-protein-reaction relationships
    into a two-column table containing reaction in column 1 and gene in column 2.

    rxn2gpr: A list of (rxnid, gpr_string) tuples
    searchkey: A regular-expression string to use to identify genes in the GPR.

    return value: a list of (reactionID, geneID) pairs. Reaction IDs with no gene IDs
    will not be included in the results.

    Default key is "fig\|\d+\.\d+\.peg\.\d+" - the expected format
    for genes in the database.'''

    genefinder = re.compile(searchkey)
    rxn2genes = []
    for tup in rxn2gpr:
        rxnid = tup[0]
        gprstring = tup[1]
        genelist = genefinder.findall(gprstring)
        for gene in genelist:
            rxn2genes.append( (rxnid, gene) )
    return rxn2genes
