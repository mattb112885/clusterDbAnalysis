#!/usr/bin/python

# This is a pipe command.
#
# Plot a scatterplot using the piped in numbers and print the result to file.
#
# 

from matplotlib import pyplot, rcParams, axis
import optparse
import fileinput
import numpy


fig = pyplot.figure(1)
fig.clf()

x = []
y = []
for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    # Skip over zeros.
    if float(spl[2]) == 0:
        continue
    x.append((float(spl[0]) + float(spl[1]))/2.0)
    y.append(float(spl[2]))

ax = fig.add_subplot(1,1,1)
ax.plot(x,y, "bo")
ax.set_xlabel("Number of elements in cluster (N=19 organisms)", {"fontsize":16})
ax.set_ylabel("Number of clusters", {"fontsize":16})
ax.set_title("Cluster size distribution", {"fontsize":16})
ax.grid(linestyle='--', b=True, which='major')

# Default number of tickmarks is too small...
ax.set_xticks(numpy.arange(0, max(x), 25))

pyplot.savefig("TEST.png", format="png")
