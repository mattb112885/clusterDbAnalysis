#!/usr/bin/python

# This is a pipe command.
# Pipe in a newick file
#
# The input argument is the "geneinfo" table containing gene IDs in the first column
# and what to replace them with in the second column.
#
# The output is a new newick file with 
# the gene names replaced (after sanitizing)

import sys
import fileinput

if not len(sys.argv) == 2:
    print "Usage: replaceGeneInNewick [Replacement table]"
    exit(2)

newickString = ''.join( [ line.strip() for line in fileinput.input("-") ] )

geneToAnnotation = {}
for line in open(sys.argv[1], "r"):
    spl = line.strip().split("\t")
    geneToAnnotation[spl[0]] = spl[1]

# Remove special characters that can confound the newick parser from the annotation strings.
# Replace with underscores
charToRemove = " \t():,'\""
for gene in geneToAnnotation:
    for char in charToRemove:
        geneToAnnotation[gene] = geneToAnnotation[gene].replace(char, "_")

# Actually replace the annotations now. We need to keep the original ID as well to make sure the names stay unique
# (grumble grumble)
for gene in geneToAnnotation:
    newickString = newickString.replace(gene, gene + "_" + geneToAnnotation[gene])

print newickString
    

