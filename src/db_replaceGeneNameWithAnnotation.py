#!/usr/bin/python

import fileinput, optparse, sqlite3, sys, re
from sanitizeString import *
from locateDatabase import *

usage="%prog [options] < infile > outfile"
description="""Look for things that look like gene IDs (fig|#.#.peg.#) in the input file
and replace them with their annotations (sanitized appropriately - so this will work in e.g. a Newick file)."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-o", "--organism", help="Include organism as part of the annotation (D: False)", action="store_true", dest="org", default=False)
parser.add_option("-k", "--keepgene", help="Keep (sanitized) gene ID as part of the annotation (D: False)", action="store_true", dest="gene", default=False)
(options, args) = parser.parse_args()

geneFinder = re.compile("fig\|\d+\.\d+\.peg\.\d+")

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

orgst = ""
genest = ""
annotest = "annotation"
if options.org:
    orgst="organism,"
if options.gene:
    genest="geneid,"

# Organism first, then gene, and finally annotation
query = "SELECT %s%s%s FROM processed WHERE processed.geneid=?;" %(orgst, genest, annotest)


for line in fileinput.input("-"):
    st = line.strip("\r\n")
    replist = geneFinder.findall(st)
    for rep in replist:
        cur.execute(query, (rep, ))

        # Get string with which to replace
        annotestr = ""
        for c in cur:
            annotestr = "_".join( [str(s) for s in c] )

        if annotestr == "":
            sys.stderr.write("WARNING: Gene id %s not found in the database - skipping...\n" %(rep))
            continue
        else:
            # Sanitize the annotation and replace it...
            annotestr = sanitizeString(annotestr, False)
            st = st.replace(rep, annotestr)
    print st

con.close()
