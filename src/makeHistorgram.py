#!/usr/bin/env python

# This is a pipe command.
# Pipe in a tab-delimited file with
# numeric entries in one of the columns
#
# Prints to stdout the beginning, end, and number of elements in the first, second, and third columns,
# respectively.
#
# Options
# -c [column number]
# -s [starting point for histogram]: Default = minimum value
# -e [ending point for histogram]: Default = maximum value
# -n Number of bins [default = 10]

from __future__ import print_function
import numpy
import optparse
import fileinput

usage = "%prog [options] < numeric_values > histogram_numbers"
description = "Calculates and prints as text the number of elements in each bin with parameters specified (does not actually make a graph for you)"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-c", "--column", help="Column number (start from 1) for numeric value (default = 1)", action="store", type="int", dest="column", default=1)
parser.add_option("-s", "--start", help="Beginning of first bin to calculate histogram for (default = min[numbers])", action="store", type="float", dest="start", default=None)
parser.add_option("-e", "--end", help="End of last bin to calculate histogram for (default = max[numbers]", action="store", type="float", dest="end", default=None)
parser.add_option("-n", "--number", help="Number of bins (default = 10)", action="store", type="int", dest="number", default=10)
(options, args) = parser.parse_args()

c = options.column - 1 # Convert to Pythonic indexes  

numarr = []
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    numarr.append(float(spl[c]))

if options.start == None:
    start = min(numarr)
else:
    start = options.start

if options.end == None:
    end = max(numarr)
else:
    end = options.end

# Calculate histogram
hist, bin_edges = numpy.histogram(numarr, bins=options.number, range=(start, end))

# bin_edges has n+1 elements - hist[i] is number of elements between bin_edges[i] and bin_edges[i+1]
for ii in range(len(hist)):
    print("%1.4f\t%1.4f\t%d" % (bin_edges[ii], bin_edges[ii+1], hist[ii]))
