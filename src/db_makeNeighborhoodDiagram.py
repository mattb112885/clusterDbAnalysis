#!/usr/bin/env python

from __future__ import print_function
import fileinput
import optparse
import os
import sqlite3
import sys

from BioPythonGraphics import *
from FileLocator import *
from sanitizeString import *

usage = "%prog runid [options] < gene_ids"

description = """Saves neighborhood diagrams for genes in a specified directory. The run ID is used to color code the neighborhood diagram. 
The color code is not necessarily consistent across different input genes but is internally consistent for neighborhoods of a given gene."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-d", "--directory", help="Directory in which to save neighborhood diagrams. Default is 'geneNeighborhoods'", action="store", 
                  dest="directory", type="str", default='geneNeighborhoods')
parser.add_option("-l", "--labeltype", help="Type of label to use. Valid types are 'aliases' or 'clusterid' (D: aliases)", action="store", dest="labeltype", type="str", default="aliases")
parser.add_option("-g", "--genecol", help="Column number for gene IDs starting from 1 (D: 1)", action="store", dest="gc", type="int", default=1)
(options, args) = parser.parse_args()

gc = options.gc - 1

if len(args) < 1:
    sys.stderr.write("ERROR: Run ID is required argument.\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

runid = args[0]

if not os.path.exists(options.directory):
    os.makedirs(options.directory)

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    geneid = spl[gc]
    diagram = makeSingleGeneNeighborhoodDiagram(geneid, runid, cur, labeltype = options.labeltype, imgfileloc = os.path.join(options.directory, sanitizeString(geneid, False)))
    sys.stderr.write("Saved result to %s\n" %(diagram))

cur.close()
