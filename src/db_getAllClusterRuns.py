#!/usr/bin/env python

# Takes no inputs.
#
# Output: A list of all run IDs from the database (if one does not desire to filter by run ID, and instead
# wishes to generate lists for ALL of the run IDs and process them later,
# one should call this function)
#
# Must be run from the root directory

import sqlite3, optparse
from ClusterFuncs import *
from FileLocator import *

header = "runid"

usage = """%prog > run_id_list

Output: """ + " ".join(header)
description = "Return list of all run IDs from the database"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

runs = getAllClusterRuns(cur)

if options.header:
    print header
print "\n".join(runs)

con.close()
