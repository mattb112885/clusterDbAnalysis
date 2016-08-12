#!/usr/bin/env python

from __future__ import print_function
import optparse

from FileLocator import *

usage = "%prog"
description = "Get the current root directory of the TIEP repository"
parser = optparse.OptionParser(usage=usage, description=description)
parser.parse_args()

print(locateRootDirectory())
