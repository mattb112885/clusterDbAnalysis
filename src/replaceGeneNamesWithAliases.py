#!/usr/bin/env python

# This is a pipe command.
# Pipe in a file in which you want the genenames replaced
#
# The input argument is the "geneinfo" table containing gene IDs in the first column
# and what to replace them with in the second column.
#
# The output is a new newick file with 
# the gene names replaced (after sanitizing)

import re
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
parser.add_option("-d", "--delimiters", help="Regex that matches anything that separates your gene IDs from other things in the string. Default is whitespace, semicolons and commas.",
                  action="store", dest="delimiters", default="[\s,;]")
(options, args) = parser.parse_args()

transfile = args[0]

# Read the translation table
geneToAnnotation = {}
for line in open(transfile, "r"):
    spl = line.strip('\r\n').split("\t")
    geneToAnnotation[spl[0]] = spl[1]

if options.nooriginal and not len(geneToAnnotation) == len(set(geneToAnnotation.values())):
    sys.stderr.write("ERROR: Replaced aliases must be unique when original names are not kept\n")
    exit(2)

# Remove special characters that can confound a Newick parser from the replacement values.
# Replace with underscores.
if not options.nosanitize:
    charToRemove = " \t():,'\""
    for gene in geneToAnnotation:
        for char in charToRemove:
            geneToAnnotation[gene] = geneToAnnotation[gene].replace(char, "_")


comp = re.compile(options.delimiters)

for originalString in fileinput.input("-"):
    originalString = originalString.strip('\r\n')
    delimiters = comp.finditer(originalString)
    delims = []
    for deli in delimiters:
        delims.append(originalString[deli.start(0):deli.end(0)])

    newstring = comp.split(originalString)
    for ii in range(len(newstring)):
        if newstring[ii] in geneToAnnotation:
            if options.nooriginal:
                newstring[ii] = geneToAnnotation[newstring[ii]]
            else:
                newstring[ii] += "_" + geneToAnnotation[newstring[ii]]

    print newstring, len(newstring)
    print delims, len(delims)
    finalstring = ''
    for ii in range(len(newstring)):
        finalstring += newstring[ii]
        if ii != len(newstring) - 1:
            finalstring += delims[ii]
            
    print finalstring
