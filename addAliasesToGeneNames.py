#!/usr/bin/python

# This is a pipe command
#
# Pipe in a data file containing annotations
# and provide the aliases file
# ( PEG ID | Alias)
#
#
# Will add aliases to the gene name (separated by underscores)
# in order specified in the aliases file.
#
# By default expects a RAW file (gene 

import optparse
import fileinput

usage = "%prog [options] aliasFile"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-a", "--annotecolumn", help="Column number for annotation (start from 1)", action="store", type="int", dest="annotec", default=8)
parser.add_option("-g", "--geneidcolumn", help="Column number for gene id (start from 1)", action="store", type="int", dest="genec", default=2)
(options, args) = parser.parse_args()

if len(args) == 0:
    print "ERROR: in %prog: Alias file must be provided"
    exit(2)

annoteColumn = options.annotec - 1
geneColumn = options.genec - 1

geneToAliases = {}
for line in open(args[0], "r"):
    # I specify '\n' because otherwise RNAs and such things as that get truncated to 12 lines. BAD.
    spl = line.strip('\r\n').split("\t")
    if spl[0] in geneToAliases:
        geneToAliases[spl[0]] += "_%s" %(spl[1])
    else:
        geneToAliases[spl[0]] = "_%s" %(spl[1])

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    if spl[geneColumn] in geneToAliases:
        spl[annoteColumn] += geneToAliases[spl[geneColumn]]
    print "\t".join(spl)
