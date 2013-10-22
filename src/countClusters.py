#!/usr/bin/env python
#
# This is a pipe command intended to take in the original cluster files from MCL
#
# Outputs the number of elements in each line (which is the number of elements in a single cluster)
#
# This can be piped into makeHistogram to get a histogram of the number of elements

import fileinput
import optparse
usage = "%prog < cluster_file > counts"
description = "Generate counts of the number of elements of each cluster from the MCL output file that has all the elements of one cluster in each row"
parser = optparse.OptionParser(usage=usage, description=description)
parser.parse_args()

for line in fileinput.input("-"):
    spl = line.split("\t")
    print len(spl)
