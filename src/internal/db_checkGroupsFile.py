#!/usr/bin/env python

from ClusterFuncs import *
from FileLocator import *
import optparse
import sqlite3
import sys

usage = "%prog"
description = "Check whether all of the organisms in the groups file are in the database."
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

groupsfile = locateGroupsFile()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

STATUS=0
for line in open(groupsfile, "r"):
    spl = line.strip("\r\n").split("\t")
    groupname = spl[0]
    organisms = spl[1].split(";")
    for org in organisms:
        try:
            organismNameToId(org, cur, issanitized = False)
        except ValueError:
            sys.stderr.write("Organism with name \"%s\" not found in database. Was this a typo? (Note - organism names with semicolons in them will cause problemd)\n" %(org))
            STATUS=1

exit(STATUS)

