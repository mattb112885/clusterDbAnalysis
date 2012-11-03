#!/usr/bin/python

import optparse
import os
import sys

usage = "%prog -o orgfile -g groupsfile -n groupname [match1] [match2] ..."
description = """Add all organims from orgfile that match [match1] [match2] ... to the groups file. 
The groups file is of the form [groupname] matchingorg1;matchingorg2;...
This function makes a group with ALL the organisms in the orgfile if nothing is specified.
"""
parser = optparse.OptionParser(usage=usage, description=description)

# FIXME - Once we have the metadata.py in place I should use it to enforce the locations of these files.
parser.add_option("-o", "--orgfile", help="organisms file (required)", action="store", type="str", dest="orgfile", default=None)
parser.add_option("-g", "--groupfile", help="groups file (required)", action="store", type="str", dest="groupfile", default=None)
parser.add_option("-n", "--groupname", help="Group name (required)", action="store", type="str", dest="groupname", default=None)
(options, args) = parser.parse_args()

if options.orgfile is None:
    sys.stderr.write("ERROR: organisms file (-o) is a required argument\n")
    exit(2)
if not os.path.exists(options.orgfile):
    sys.stderr.write("ERROR: Specified organisms file %s does not exist\n" %(options.orgfile))
if options.groupfile is None:
    sys.stderr.write("ERROR: groups file (-g) is a required argument\n")
    exit(2)
if options.groupname is None:
    sys.stderr.write("ERROR: group name (-n) is a required argument\n")
    exit(2)
if len(args) == 0:
    sys.stderr.write("ERROR: At least one key to match must be provided\n")
    exit(2)

orglist = set()
for line in open(options.orgfile, "r"):
    spl = line.strip("\r\n").split("\t")
    orglist.add(spl[0])

group2orglist = {}
if os.path.exists(options.groupfile):
    for line in open(options.groupfile, "r"):
        spl = line.strip("\r\n").split("\t")
        # I use set so that order of organisms doesn't matter...
        group2orglist[spl[0]] = set(spl[1].split(";"))

if options.groupname in group2orglist:
    sys.stderr.write("ERROR: Specified group name %s already exists in groups file \n" %(options.groupname))
    exit(2)

# Look for matching organisms (enforcing uniqueness in case we get multiple possible matches that all match the same organism)
matchingorgs = set()
if len(args) == 0:
    sys.stderr.write("No matches provided - will use all the organisms in the named group...\n")
    matchingorgs = orglist
else:
    for arg in args:
        for org in orglist:
            if arg.lower() in org.lower():
                matchingorgs.add(org)

if len(matchingorgs) == 0:
    sys.stderr.write("ERROR: No organisms matched the specified queries.\n")
    exit(2)

# Don't cluster the same group with multiple names.
if matchingorgs in group2orglist.values():
    sys.stderr.write("ERROR: The list of organisms matching your query is already present in the groups file under a different name\n")
    exit(2)

# Now that our sanity checks have passed...
linetoadd = "%s\t%s\n" %(options.groupname, ";".join(sorted(list(matchingorgs))))
groupfid = open(options.groupfile, "a+")
groupfid.write(linetoadd)
groupfid.close()
sys.stderr.write("New line written to groups file\n")
