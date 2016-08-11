#!/usr/bin/env python

import fileinput, math, optparse, sys

okmethods = ["pctid", "logevalue", "minbit", "maxbit"]

usage = "%prog [options] < BLAST_results > distance_matrix"
description = """Turn a table of BLAST results into a distance matrix.
The distance matrix will be suitable for conversion into a heatmap with
plotHeatmap.py
It is suggested that the user use this with the results of
db_getBlastResultsBetweenSpecificGenes.py so that there are blast results
available for ALL pairs of genes in a specified list."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--method", help="Scoring metric to use (D=maxbit). OK methods: %s" %(" ".join(okmethods)), action="store", type="str", dest="method", default="maxbit")
parser.add_option("-d", "--default", help="Default score to use in case of missing pairs. (D: None)", action="store", type="float", dest="default", default=None)
(options, args) = parser.parse_args()

if not options.method in okmethods:
    sys.stderr.write("ERROR: Provided method %s is not in list of recognized scoring methods: %s\n" %(options.method, " ".join(okmethods) ) )
    exit(2)

# Keep unique until the end
queries = set()
targets = set()
# Dict from (query,target) pairs to the score between them
qt_to_score = dict()

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    if options.method == "pctid":
        score = float(spl[3])
    elif options.method == "logevalue":
        score = math.log10(float(spl[10]))
    elif options.method == "minbit":
        score = float(spl[11])/min(float(spl[12]), float(spl[13]))
    elif options.method == "maxbit":
        score = float(spl[11])/min(float(spl[12]), float(spl[13]))
    queries.add(spl[0])
    targets.add(spl[1])
    if (spl[0], spl[1]) in qt_to_score:
        sys.stderr.write("WARNING: Multiple instances of the pair (%s, %s) in your list - only the first will be kept\n" %(spl[0], spl[1]))
        continue
    qt_to_score[ (spl[0], spl[1]) ] = score

queries = sorted(list(queries))
targets = sorted(list(targets))

print("\t" + "\t".join(targets))

for query in queries:
    row = query
    for target in targets:
        if (query, target) in qt_to_score:
            row = row + "\t" + str(qt_to_score[(query, target)])
        else:
            row = row + "\t" + str(options.default)
    print(row)
