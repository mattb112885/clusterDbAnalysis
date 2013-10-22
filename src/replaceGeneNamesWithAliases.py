#!/usr/bin/env python

# This is a pipe command.
# Pipe in a file in which you want the genenames replaced
#
# The input argument is the "geneinfo" table containing gene IDs in the first column
# and what to replace them with in the second column.
#
# The output is a new newick file with 
# the gene names replaced (after sanitizing)

import sys
import fileinput
import optparse

description="Replace one name of a gene with another name as given in a specified replacement table"
usage="%prog (options) replacement_table < file > file_with_geneids_replaced"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--nooriginal", help="Set this flag to NOT keep the original name. To use this the new aliases must be unique (D=False)", 
                  action="store_true", dest="nooriginal", default=False)
parser.add_option("-s", "--nosanitize", help="""Set this flag to NOT sanitize the replacement names. 
Sanitized characters do not include pipes or periods so that this function can be used to put gene IDs in without this function""", 
                  action="store", dest="nosanitize", default=False)
(options, args) = parser.parse_args()

transfile = args[0]

newickString = '\n'.join( [ line.strip('\r\n') for line in fileinput.input("-") ] )
geneToAnnotation = {}
for line in open(transfile, "r"):
    spl = line.strip('\r\n').split("\t")
    geneToAnnotation[spl[0]] = spl[1]

if options.nooriginal and not len(geneToAnnotation) == len(set(geneToAnnotation.values())):
    sys.stderr.write("ERROR: Replaced aliases must be unique when original names are not kept\n")
    exit(2)

# Remove special characters that can confound the newick parser from the annotation strings.
# Replace with underscores
if not options.nosanitize:
    charToRemove = " \t():,'\""
    for gene in geneToAnnotation:
        for char in charToRemove:
            geneToAnnotation[gene] = geneToAnnotation[gene].replace(char, "_")

# Actually replace the annotations now. We need to keep the original ID as well to make sure the names stay unique
# (grumble grumble)
for gene in geneToAnnotation:
    if options.nooriginal:
        newickString = newickString.replace(gene, geneToAnnotation[gene])
    else:
        newickString = newickString.replace(gene, gene + "_" + geneToAnnotation[gene])

print newickString
    

