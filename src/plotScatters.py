#!/usr/bin/python

# This is a pipe command.
#
# Plot a scatterplot using the piped in numbers and print the result to file.
# Plots using relatively nice font sizes, dot sizes, etc...
# 

from matplotlib import pyplot, rcParams, axis
import optparse
import fileinput
import numpy
import math
import sys

# Estimate bounds that will make for an aesthetically pleasing plot.
# As a heuristic, I say the start and end have to be divisible by 2
# (within the order of magnitude)
def estimateBounds(arr):
    # Get order of magnitude estimate for each datum
    oom = []
    for a in arr:
        if a == 0:
            oom.append(0)
        else:
            oom.append(int(math.log10(abs(a))))

    # Warn if there are too many order of magnitude differences - should maybe be using a log plot!
    if max(oom) - min(oom) > 2:
        sys.stderr.write("WARNING: Wide range of orders of magnitude detected. If the plot does not look good, try a log plot!\n")

    # Standardize the order of magnitude of the bounding points
    minval = min(arr)
    maxval = max(arr)
    mintrans = minval * 10**-max(oom)
    maxtrans = maxval * 10**-max(oom)

    if maxtrans > 10 or mintrans < -10:
        sys.stderr.write("INTERNAL ERROR: I'm a retard\n")
        exit(2)

    # It is generally acceptable to have bounds in multiples of 2, 5, or 10.
    okbounds = [ (0,2), (0,5), (0,10), (-10, 0), (-5, 0), (-2,0), (-10,10), (-5,5), (-2,2) ]
    mindev = 100
    bestbounds = None
    for bounds in okbounds:
        # Minimum value of x must be bigger than the minimum on the plot
        if bounds[0] > mintrans:
            continue
        # Maximum value of x must be less than the maximum on the plot
        if bounds[1] < maxtrans:
            continue
        diff = (mintrans - bounds[0]) + (bounds[1] - maxtrans)
        if diff < mindev:
            mindev = diff
            bestbounds = bounds

    # Re-transform to the original order of magnitude and return the result.
    return bestbounds[0]*10**max(oom), bestbounds[1]*10**max(oom)

usage="%prog [options] < x,y"
description="User interface to generate a single scatterplot. The default settings should make a pretty nice plot but I provide easy-to-understand options to change different aspects of the plot"
parser = optparse.OptionParser(usage=usage, description=description)
# Input Data options
parser.add_option("--xcol", help="Column number for x, starting from 1 (D = 1)", action="store", type="int", dest="xcol", default=1)
parser.add_option("--ycol", help="Column number for y, starting from 1 (D = 2)", action="store", type="int", dest="ycol", default=2)
# Output data options
parser.add_option("--noshow", help="Do not display plot (displaying requires X server)", action="store_false", dest="show", default=True)
parser.add_option("--png", help="Save plot as the specified file in png format", action="store", type="str", dest="pngout", default=None)
# Plot options
parser.add_option("--minx", help="Minimum value for x-axis (Default: Auto-estimate)", action="store", type="float", dest="minx", default=None)
parser.add_option("--maxx", help="Maximum value for x-axis (Default: Auto-estimate)", action="store", type="float", dest="maxx", default=None)
parser.add_option("--miny", help="Maximum value for y-axis (Default: Auto-estimate)", action="store", type="float", dest="miny", default=None)
parser.add_option("--maxy", help="Maximum value for y-axis (Default: Auto-estimate)", action="store", type="float", dest="maxy", default=None)
parser.add_option("--fontsize", help="Font size for axis labels (Default: 16)", action="store", type="int", dest="fontsize", default=16)
parser.add_option("--xlog", help="Put x-axis on log scale (Default: Linear scale)", action="store_true", dest="xlog", default=False)
parser.add_option("--ylog", help="Put y-axis on log scale (Default: Linear scale)", action="store_true", dest="ylog", default=False)
parser.add_option("--xlabel", help="Label for X-axis (default: warning - no label", action="store", type="str", dest="xlabel", default="WARNING: No X label provided")
parser.add_option("--ylabel", help="Label for Y-axis (default: warning - no label", action="store", type="str", dest="ylabel", default="WARNING: No Y label provided")
parser.add_option("--title", help="Label for title (default: warning - no title", action="store", type="str", dest="title", default="WARNING: No title provided")
parser.add_option("--connect", help="If set, connect the dots, otherwise they're just sepatate", action="store_true", dest="connect", default=False)
parser.add_option("--connectdots", help="If set, connect dots and display the dots in addition to the connecting line (implies --connect)", action="store_true", dest="connectdots", default=False)
parser.add_option("--color", help="Set color of dots (default: as chosen by matplotlib)", action="store", type="str", dest="color", default="blue")
(options, args) = parser.parse_args()

# If you specify one bound you have to specify both
if (options.minx == None and options.maxx != None) or (options.minx != None and options.maxx == None):
    sys.stderr.write("Only one bound specified on X. If you specify one bound on X, you must specify both\n")
    exit(2)
if (options.miny == None and options.maxy != None) or (options.miny != None and options.maxy == None):
    sys.stderr.write("Only one bound specified on Y. If you specify one bound on Y, you must specify both\n")
    exit(2)

xc = options.xcol - 1
yc = options.ycol - 1

xyp = []
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    xyp.append( (float(spl[xc]), float(spl[yc])))

# Sort by x to avoid all kinds of weird stuff when connecting dots
xyp = sorted(xyp)

if options.connect or options.connectdots:
    dispfmt = "-"
else:
    dispfmt = "o"

x = [ s[0] for s in xyp ]
y = [ s[1] for s in xyp ]

pyplot.plot(x,y, dispfmt, color=options.color)
if options.connectdots:
    pyplot.plot(x,y,"o", color=options.color)

pyplot.xlabel(options.xlabel, fontsize=options.fontsize)
pyplot.ylabel(options.ylabel, fontsize=options.fontsize)
pyplot.title(options.title, fontsize=options.fontsize)

if options.xlog:
    if min(x) <= 0:
        sys.stderr.write("ERROR: Specified X-data contains negative values so cannot use log plot on X axis\n")
        exit(2)
    pyplot.xscale('log')
if options.ylog:
    if min(y) <= 0:
        sys.stderr.write("ERROR: Specified Y-data contains negative values so cannot use log plot on Y axis\n")
        exit(2)
    pyplot.yscale('log')

# These have to be done last because lots of other options seem to clobber it (including making log scale on the Y axis clobbering our changes to the X axis)
if options.minx == None and options.maxx == None and not options.xlog:
    xbounds = estimateBounds(x)
    pyplot.xlim(xbounds)
elif not options.xlog:
    pyplot.xlim( (options.minx, options.maxx) )

if options.miny == None and options.maxy == None and not options.ylog:
    ybounds = estimateBounds(y)
    pyplot.ylim(ybounds)
elif not options.xlog:
    pyplot.ylim( (options.miny, options.maxy) )
    
if options.pngout != None:
    pyplot.savefig(options.pngout, format="png")

if options.show:
    pyplot.show()
