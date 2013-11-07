#!/usr/bin/env python

# Given a list of organisms to stdin and a cluster run,
# identify clusters within that run that only have representatives
# within the specified list of organisms

import optparse
import sys
from CoreGeneFunctions import *

usage="%prog (-a|-n|-p|-s|-y) [options] run_id < organism_list > cluster_run_id_list"
description="""Find clusters with a paritcular quality relative to the list of organisms you specified.
Note: To find conserved gene clusters use -a
To find core gene clusters for a particular group, use both -a and -u
To find core genes only in a parituclar group (to the exclusion of all the others in that cluster run), use -a, -u, and -s
To find clusers that exclude all the specified organisms use -n
Using only -u or various contradictory combinations will result in an error.
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-a", "--all", help="Only include clusters that have at least one representative in ALL of the specified organisms", action="store_true", dest="all", default=False)
parser.add_option("-y", "--any", help="Only include clusters that have at least one representative in AT LEAST ONE of the specified organisms", action="store_true", dest="any", default=False)
parser.add_option("-s", "--only", help="Only include clusters that ONLY has matches in the specified organisms", action="store_true", dest="only", default=False)
parser.add_option("-n", "--none", help="Only include clusters that have DOES NOT have a representative in the specified organisms", action="store_true", dest="none", default=False)
parser.add_option("-u", "--uniq", help="Only include clusters that contain exactly ONE representative in any matching organisms (D: Any number)", action="store_true", dest="uniq", default=False)
parser.add_option("-o", "--orgcol", help="Column number for organism starting from 1 (D=1)", action="store", type="int", dest="oc", default=1)
parser.add_option("-r", "--sanitized", help="If specified, the names in the input file have been sanitized (with sanitizeString.py) (D: False)", action="store_true", dest="sanitized", default=False)
parser.add_option("-p", "--pct_cutoff", help="Percentage of organisms in the specified list of organisms that must have a member to be included (incompatible with -a, -y, and -n but can be used with -s)",
                  action="store", type="float", dest="pct_cutoff", default=None)
#parser.add_option("-m", "--minorgs", help="Minimum number of organisms in clusters to be included (D=no minimum)", action="store", type="int", dest="minorg", default=None)
(options, args) = parser.parse_args()

if not len(args) == 1:
    sys.stderr.write("ERROR: Run ID must be provided\n")
    exit(2)

if options.all and options.none:
    sys.stderr.write("ERROR: ALL and NONE options are contradictory\n")
    exit(2)

if options.any and options.none:
    sys.stderr.write("ERROR: ANY and NONE options are contradictory\n")
    exit(2)

if options.only and options.none:
    sys.stderr.write("ERROR: ONLY and NONE options are contradictory\n")
    exit(2)

if not (options.only or options.all or options.any or options.none or options.pct_cutoff is not None):
    sys.stderr.write("ERROR: At least one of -a, -y, -s, -n, -p must be specified\n")
    exit(2)

oc = options.oc - 1

orglist = []
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    orglist.append(spl[oc])

clusterrun_list = findGenesByOrganismList(orglist, args[0], 
                                          sanitized = options.sanitized, 
                                          any_org = options.any,
                                          all_org = options.all,
                                          only_org = options.only,
                                          none_org = options.none,
                                          uniq_org = options.uniq,
                                          pct_cutoff = options.pct_cutoff
                                          )

for cr in clusterrun_list:
    print "%s\t%s" %(cr[0], cr[1])
