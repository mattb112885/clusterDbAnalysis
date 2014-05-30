#!/usr/bin/env python

import fileinput, optparse, sqlite3, sys, re
from FileLocator import *
from sanitizeString import *

usage="""
%prog [options] < infile > outfile
Replace with organism only: use -o
Replace with annotation only: use -a
Replace with organism and annotation and keep original (sanitized) gene id: use -a -o -k
"""

description="""Look for things that look like gene IDs (fig|#.#.peg.#) in the input file
and replace them with annotation and\or organism name, properly sanitized to work in a Newick file.
Also works if the input has sanitized gene IDs (fig_#_#_peg_#). """

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-a", "--annote", help="Include annotation (D: False)", action="store_true", dest="ann", default=False)
parser.add_option("-o", "--organism", help="Include organism as part of the annotation (D: False)", action="store_true", dest="org", default=False)
parser.add_option("-k", "--keepgene", help="Keep (sanitized) gene ID as part of the annotation (D: False)", action="store_true", dest="gene", default=False)
(options, args) = parser.parse_args()

if not (options.ann or options.org or options.gene):
    sys.stderr.write("ERROR: Must specify -a, -o, or -k in db_replaceNameWithAnnotation\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

wanted = []
if options.org:
    wanted.append("organism")
if options.gene:
    wanted.append("geneid")
if options.ann:
    wanted.append("annotation")
selstr = ",".join(wanted)
query = "SELECT %s FROM processed WHERE processed.geneid=?;" %(selstr)

geneFinder = re.compile("fig\|\d+\.\d+\.peg\.\d+")

for line in fileinput.input("-"):

    st = line.strip("\r\n")
    st = unsanitizeGeneId(st)
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
