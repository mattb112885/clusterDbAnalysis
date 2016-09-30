#!/usr/bin/env python

from __future__ import print_function
import optparse
import sys

usage = "%prog aliases file file_to_replace"
description = """Given a file containing aliased gene names,
replaces those aliases with the corresponding gene ID (e.g. fig|...).
I used this to take GPRs from a model (with corresponding locus tags)
and translate them into a form usable in the database.
WARNING: This is a dumb script; anything that shouldnt be an alias
but that matches an alias (like 'for') will be replaced. I highly suggest
using this only for locus tags and making sure no other aliases match
the locus tags..."""
parser = optparse.OptionParser(usage=usage, description=description)
(opts, args) = parser.parse_args()

if len(args) < 2:
    sys.stderr.write("ERROR: Both aliases file and the file with aliases to replace are required arguments to replaceAliasesWithGeneNames.py\n")
    exit(2)

alias2gene = {}
for line in open(args[0]):
    spl = line.strip("\r\n").split("\t")
    alias2gene[spl[1]] = spl[0]

for line in open(args[1]):
    spl = line.strip("\r\n").split("\t")
    repl = spl[1]
    for alias in alias2gene:
        if alias in repl:
            repl = repl.replace(alias, alias2gene[alias])
    print("%s\t%s" %(spl[0], repl))
