#!/usr/bin/python

from __future__ import print_function
import optparse
import os
import re
import shutil
import sys
from FileLocator import *

usage = "%prog (-f|-l) -r runid [file]"
description = '''Import a clustering result obtained outside ITEP into the ITEP database 

Given a file containing a clustering of genes with ITEP IDs,
import the clustering into ITEP. Clustering format can be either in "line" format (-l) or
"flat" format (-f).

"Line" format is a tab-deimited file of gene IDs with one line per cluster. In the following example,
gene IDs 1-3 would be in one cluster and gene IDs 4 and 5 would be in a different cluster:

[geneid1] [geneid2] [geneid3]
[geneid4] [geneid5]

In "flat" format, there are three tab-delimited columns: a run ID (specifying how the cluster was created),
a cluster ID (typically an integer) and a gene ID. The following example represents the equivalent clustering
to the example above, in flat format:

[runid] 1 [geneid1]
[runid] 1 [geneid2]
[runid] 1 [geneid3]
[runid] 2 [geneid4]
[runid] 2 [geneid5]

There should be no headers in either type of file.

If MCL format is provided, the run ID is either what is user-specified or (if not specified) is assigned to be
the same name as the input file. Cluster IDs are generated as sequential integerrs (same as what is done when 
running setup_step2.sh).

If flat format is provided, you should have already assigned a run ID (or IDs) and therefore you will get an
error if attempting to assign a different run ID.

Gene IDs must match IDs in the ITEP database (e.g. fig|\d+\.\d+\.peg\.\d+). You will get foreign key exceptions if
this is not the case and the clusters will not be loaded.
'''

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--flatfile", help="Specifies that the file is in flat format (see description). Either -f or -l must be specified",
                  action="store_true", dest="flatfile", default=False)
parser.add_option("-l", "--linefile", help="Specifies that the file is in line format (see description). Either -f or -l must be specified",
                  action="store_true", dest="linefile", default=False)
parser.add_option("-r", "--runid", help="Specify a run ID to use when processing the input file. Can only be used with -l. By default, we use the name of the input file as the run ID.",
                  action="store", dest="runid", default=None)
(options, args) = parser.parse_args()

if len(args) < 1:
    sys.stderr.write("ERROR: Input file must be specified (see -h for usage instructions)\n")
    exit(2)

if options.flatfile and options.linefile:
    sys.stderr.write("ERROR: Cannot specify both flat (-f) and line (-l) types of input files.\n")
    exit(2)

if not(options.flatfile or options.linefile):
    sys.stderr.write("ERROR: Must specify either flat (-f) or line (-l) types of input files.\n")
    exit(2)

if options.flatfile and options.runid is not None:
    sys.stderr.write("ERROR: Cannot specify a new run id when using a flat file because run ID should already be in the first column of the flat file.\n")
    exit(2)

infile = args[0]

if not os.path.isfile(infile):
    sys.stderr.write("ERROR: Input file %s does not exist!\n" %(infile))
    exit(2)

# If cluster and flatcluster directories do not exist, create them
# (this is necessary if setup_step2.sh hasn't been run yet, i.e. if a user ONLY wants
# to use clusters generated externally, and not use MCL from ITEP's scripts.
rootdir = locateRootDirectory()
clusterdir = os.path.join(rootdir, "clusters")
flatclusterdir = os.path.join(rootdir, "flatclusters")

if not os.path.exists(clusterdir):
    os.path.mkdir(clusterdir)
if not os.path.exists(flatclusterdir):
    os.path.mkdir(flatclusterdir)

# Check the inpurt formats. Also, 
# Get the run ID. How we get this depends on the type of file:
# linefile - get from filename (default) or from options
# flatfile - get from the file
# We will create copy to a file with this name but will throw an error if there already
# is a file of that name there.
genefinder = re.compile("^fig\|\d+\.\d+\.peg\.\d+$")
if options.flatfile:
    with open(infile, 'r') as f:
        first_line = f.readline()
        spl = first_line.strip("\r\n").split("\t")
        if len(spl) != 3:
            sys.stderr.write("ERROR: -f specified but file is not in expected format (did you mean to specify -l?). Expected format has three columns: run ID, cluster ID and gene ID (see -h)\n")
            exit(2)
        runid = spl[0]
        geneid = spl[2]
        if genefinder.match(geneid) is None:
            sys.stderr.write("ERROR: Gene IDs in the input file were not in expected format (must match IDs in the ITEP database). Offending gene ID: %s\n" %(geneid))
            exit(2)
        if genefinder.match(runid) is not None:
            sys.stderr.write("ERROR: The column expected to be run ID matched gene ID - it is very likely that you meant to specify -f instead or that the columns are in the wrong order\n")
            exit(2)
        destdir = os.path.join(rootdir, "flatclusters", runid)
else:
    with open(infile, 'r') as f:
        first_line = f.readline()
        spl = first_line.strip("\r\n").split("\t")
        geneid = spl[0]
        if genefinder.match(geneid) is None:
            sys.stderr.write("ERROR: Gene IDs in the input file were not in expected format (must match IDs in the ITEP database). Offending gene ID: %s\n" %(geneid))
            sys.stderr.write("(did you mean to specify -f?)\n")
            exit(2)
    if options.runid is None:
        runid = infile
    else:
        runid = options.runid
    destdir = os.path.join(rootdir, "clusters", runid)

if os.path.isfile(destdir):
    sys.stderr.write("ERROR: Destination cluster file %s already exists - please delete the file and rerun if you really want to replace it.\n" %(destdir))
    exit(2)
else:
    shutil.copyfile(infile, destdir)

# Run the import script which flattens cluster files, concatenates them and runs the SQL to import it into the database.
os.chdir(rootdir)
os.system("./importAllClusters.sh")
