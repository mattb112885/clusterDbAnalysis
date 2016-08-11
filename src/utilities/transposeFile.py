#!/usr/bin/env python

import fileinput
import numpy
import optparse

usage = "%prog < delimited_file > transposed_delimited_file"
description = """Transpose all the stuff in a file (I wrote this to get around Excel's shockingly small column limit)"""

parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

X = numpy.array( [ tuple(line.strip("\r\n").split("\t")) for line in fileinput.input("-") ] )

X = X.transpose()

for ii in X:
    print("\t".join(ii))
