#!/usr/bin/python
#
# This is a pipe command intended to take in the original cluster files from MCL
#
# Outputs the number of elements in each line (which is the number of elements in a single cluster)
#
# This can be piped into makeHistogram to get a histogram of the number of elements

import fileinput

for line in fileinput.input("-"):
    spl = line.split("\t")
    print len(spl)
