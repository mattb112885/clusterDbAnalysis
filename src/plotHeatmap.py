#!/usr/bin/python

import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
import matplotlib.pyplot as pyplot
import matplotlib
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
description = "Generate a heat map for an input of numeric data. It comes with a rich set of options."
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--outfile", help="Name of output file (D: Just display graph on screen)", action="store", type="str", dest="outfile", default=None)
parser.add_option("-r", "--rowdendrogram", help="Make a row dendrogram (D: Heatmap only)", action="store_true", dest="rowdendrogram", default=False)
parser.add_option("-c", "--columndendrogram", help="Make a column dendrogram (D: Heatmap only)", action="store_true", dest="columndendrogram", default=False)
parser.add_option("-w", "--rowlabels", help="Specify this flag if the input file contains row labels (D: Label by position in original file)", action="store_true", dest="rowlabels", default=False)
parser.add_option("-o", "--columnlabels", help="Specify this flag if the input file contains column labels (D: Label by position in original file)", action="store_true", dest="columnlabels", default=False)
parser.add_option("-d", "--distancemetric", help="Distance metric between rows\columns of the input matrix (D: euclidean)", action="store", type="str", dest="distancemetric", default="euclidean")
parser.add_option("-m", "--clustermethod", help="Clustering method to use (D:complete - complete-linkage clustering)", action="store", type="str", dest="clustermethod", default="complete")
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

# Figure setup
fig = pyplot.figure(1)

# First dendrogram (left)
# Note - format is fraction from [ left, bottom, width, height] 
# bottom-left is (0,0)
if options.rowdendrogram:
    ax = pyplot.axes([0.05,0.1,0.08,0.6], frameon=False)
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
ax = pyplot.axes([0.3,0.1,0.6,0.6])
ax.set_xticks([])
ax.set_yticks([])
im = pyplot.pcolor(mat, cmap="gray")

# Colorbar
ax = pyplot.axes([0.91,0.1,0.02,0.8], frameon=False)
pyplot.colorbar(im, cax=ax)

if options.outfile is not None:
    pyplot.savefig(options.outfile)

pyplot.show() 
