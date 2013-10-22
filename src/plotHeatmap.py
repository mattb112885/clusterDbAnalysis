#!/usr/bin/env python

import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
import matplotlib.pyplot as pyplot
import numpy
import fileinput, optparse, sys

def displayDendrogram(mat, metric, method, labels, dimension):
    # Dimension = 1: Row-wise
    # Dimension = 2: Column-wise
    orientation = ""
    if dimension == 1:
        tmat = mat
        orientation = "right"
    elif dimension == 2:
        tmat = numpy.transpose(mat)
        orientation = "top"
    else:
        raise IOError
    # Distance matrix and clustering
    Y = dist.pdist(tmat, metric)
    lnk = hier.linkage(Y, method=method)
    dend = hier.dendrogram(lnk, orientation=orientation, labels=labels)
    # Sort (original) matrix according to clustering results
    if dimension == 1:
        mat = mat[dend["leaves"], :]
    elif dimension == 2:
        mat = mat[:, dend["leaves"] ]
    return mat

####################

usage = "%prog [options] < Input_tsv"
description = """Generate a heat map for a tab-delimited input of numeric data. It comes with a rich set of options.
If row labels or column labels are present in the file you must specify -w (row) or -o (column) or the tab-delimited file will fail to parse.
Labels are only actually put on the graph if the appropriate dendrogram is desired (-r for row or -c for column).
Distance method is passed directly to scipy.spatial.distance.pdist() - see docs on that function for details on valid parameters.
Cluster method is passed directly to scipy.cluster.hierarchy.linkage() - see docs on that function for details on valid parameters.
Color map is one of the colormaps recognized by pyplot - see http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--outfile", help="Name of output file (D: Just display graph on screen)", action="store", type="str", dest="outfile", default=None)
parser.add_option("-r", "--rowdendrogram", help="Make a row dendrogram (D: Heatmap only)", action="store_true", dest="rowdendrogram", default=False)
parser.add_option("-c", "--columndendrogram", help="Make a column dendrogram (D: Heatmap only)", action="store_true", dest="columndendrogram", default=False)
parser.add_option("-w", "--rowlabels", help="Specify this flag if the input file contains row labels (D: Label by position in original file)", action="store_true", dest="rowlabels", default=False)
parser.add_option("-o", "--columnlabels", help="Specify this flag if the input file contains column labels (D: Label by position in original file)", action="store_true", dest="columnlabels", default=False)
parser.add_option("-d", "--distancemetric", help="Distance metric between rows\columns of the input matrix (D: euclidean)", action="store", type="str", dest="distancemetric", default="euclidean")
parser.add_option("-m", "--clustermethod", help="Clustering method to use (D:complete - complete-linkage clustering)", action="store", type="str", dest="clustermethod", default="complete")
parser.add_option("-a", "--colormapscheme", help="Color map coloring scheme to use (D: gray)", action="store", type="str", dest="colormapscheme", default="gray")
parser.add_option("-x", "--minscore", help="Minimum score to display on color map (D: Scale to provided values)", action="store", type="float", dest="minscore", default=None)
parser.add_option("-y", "--maxscore", help="Maximum score to display on color map (D: Scale to provided values)", action="store", type="float", dest="maxscore", default=None)
(options, args) = parser.parse_args()

rowlabels = None
columnlabels = None

mat = []

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    if options.columnlabels and columnlabels is None:
        # Top-left corner is empty if both column and row labels are present
        if options.rowlabels:
            columnlabels = spl[1:]
        else:
            columnlabels = spl[0:]
        continue
    # Row label in first column (the continue above ensures that we don't try to get a
    # row label from the column label row)
    if options.rowlabels:
        if rowlabels is None:
            rowlabels = [ spl[0] ]
        else:
            rowlabels.append(spl[0])
        mat.append(spl[1:])
    else:
        mat.append(spl[0:])

# Need to convert matrix to float or pyplot will whine.
mat = numpy.array(mat, dtype="float32")

# Sanity checks
flatmat = numpy.ndarray.flatten(mat)
if options.minscore is not None and min(flatmat) < options.minscore:
    sys.stderr.write("WARNING: The requested minimum score on the color bar would leave off some data in the provided data file. Is this a mistake?\n")
    sys.stderr.write("Requested minimum: %1.4f ; Data minimum: %1.4f\n" %(options.minscore, min(flatmat)))
if options.maxscore is not None and max(flatmat) > options.maxscore:
    sys.stderr.write("WARNING: The requested maximum score on the color bar would leave off some data in the provided data file. Is this a mistake?.\n")
    sys.stderr.write("Requested maximum: %1.4f ; Data maximum: %1.4f\n" %(options.maxscore, max(flatmat)))
if options.rowlabels and not options.rowdendrogram:
    sys.stderr.write("WARNING: Without specifying -r (--rowdendrogram), row labels will not be printed on the graph\n")
if options.columnlabels and not options.columndendrogram:
    sys.stderr.write("WARNING: Without specifying -c (--columndendrogram), column labels will not be printed on the graph\n")

# Figure setup
fig = pyplot.figure(1, figsize=(18,13))

# First dendrogram (left)
# Note - format is fraction from [ left, bottom, width, height] 
# bottom-left is (0,0)
if options.rowdendrogram:
    ax = pyplot.axes([0.02,0.02,0.08,0.6], frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    mat = displayDendrogram(mat, options.distancemetric, options.clustermethod, rowlabels, 1)

# Second dendrogram (top)
if options.columndendrogram:
    ax = pyplot.axes([0.3,0.75,0.6,0.2], frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    mat = displayDendrogram(mat, options.distancemetric, options.clustermethod, columnlabels, 2)

# Heatmap
ax = pyplot.axes([0.3,0.02,0.6,0.6])
ax.set_xticks([])
ax.set_yticks([])
im = pyplot.pcolor(mat, cmap=options.colormapscheme, vmin=options.minscore, vmax=options.maxscore)

# Colorbar
ax = pyplot.axes([0.91,0.02,0.02,0.6], frameon=False)
pyplot.colorbar(im, cax=ax)

if options.outfile is not None:
    pyplot.savefig(options.outfile, dpi=300)

pyplot.show() 
