#!/usr/bin/env python

import fileinput
import optparse
import operator
import sqlite3
import sys
from ClusterFuncs import *
from FileLocator import *
from getSequenceRegion import *

usage = "%prog [options] < TBLASTN_output > bad_mutation_list"
description = """Given a TBLASTN output from db_TBlastN_wrapper,
attempt to automatically identify frameshifts and fragments of the same
gene that appear on separate lines.

The criteria for a frameshift is if the TBLASTN interval is on the same
strand in different frames in two differnet hits and they are within INTERVAL
base pairs of each other.

The criteria for a nonsense mutation is if a single hit has an internal
stop codon.

The criteria for an insertion is if there are two hits within INTERVAL
base pairs of each other.

All of these require manual curation to ensure that alternative causes
such as gene duplications and two distantly-related but close-proximity genes
are not the cause for the observed anomolies.

This function only works if we can access the database and pull out contigs.
It won't work for TBLASTN results from arbitrary fastn files"""

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-i", "--interval", 
                  help="Interval between two consecutive hits within which we will assume that the two hits belong to the same ancestral gene (D=1000 basepairs)",
                  action="store", type="int", dest="interval", default=1000)
parser.add_option("-n", "--nonsense_pct",
                  help = """Only report nonsense mutations if more than -n percent of the amino acids in the homnologous regions are found after the stop codon (cut off)
(D = 0: report any stop codons that are found)""",
                  action = "store", type="int", dest="nonsense_pct", default = 0)
parser.add_option("-f", "--falsepositives",
                  help = "Print very likely false positive issue identifications (multiple TBLASTN hits with exactly the same coordinates) - default is to ignore them",
                  action = "store_true", dest="falsepositives", default=False)

(options, args) = parser.parse_args()

# result_dict[contig][query_gene] = [ ( start, stop strand ), (start, stop, strand) ... ]
result_dict = {}

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    query_gene = spl[0]
    target_contig = spl[2]
    target_org = spl[3]
        
    start = int(spl[4])
    stop = int(spl[5])
    strand = int(spl[10])
    target_tuple = ( start, stop, strand )

    if target_contig in result_dict:
        if query_gene in result_dict[target_contig]:
            result_dict[target_contig][query_gene].append( target_tuple )
        else:
            result_dict[target_contig][query_gene] = [ target_tuple ]
    else:
        result_dict[target_contig] = {}
        result_dict[target_contig][query_gene] = [ target_tuple ]

contig_q = """SELECT seq FROM contigs WHERE contig_mod = ?"""

for contig in result_dict:
    # Get the contig sequence so that we can ask whether or not there is a stop codon
    # within the hit or not...
    res = cur.execute(contig_q, (contig, ))
    whole_contig_seq = None
    for c in res:
        whole_contig_seq  = c[0]
    if whole_contig_seq is None:
        sys.stderr.write("ERROR: I screwed up. contig name = %s \n" %(contig))
        exit(2)

    for query in result_dict[contig]:

        ls = result_dict[contig][query]
        # Sort by start then by stop
        # This isn't a rigorous approach - implementing with interval trees would be better so that we
        # don't miss any overlaps.
        ls = sorted(ls, key=operator.itemgetter(0,1) )

        # Test for multiple hits within DIST of each other which could indicate insertions or frameshifts
        for ii in range(len(ls) - 1):
            # Type of anomoly (will be filled in later and has no issues if stays as NONE)
            TYPE = None

            first_comp = ls[ii]
            second_comp = ls[ii+1]
            # As an artefact of how I did my overlapping gene analysis some TBLASTN hits are identical
            # and it is very likely that these are false positives.
            # By default I throw these away. If you really want to look at them you can set a flag to keep them.
            if first_comp == second_comp and not options.falsepositives:
                continue

            # Now thet we know they are actually DIFFERENT hits...
            dist = min( abs(first_comp[0] - second_comp[0]),
                        abs(first_comp[0] - second_comp[1]),
                        abs(first_comp[1] - second_comp[0]),
                        abs(first_comp[1] - second_comp[1]) 
                        )

            # Attempt to guess the type
            # These will require further curation e.g. to ensure that each half of the frameshift
            # matches a different part of the query protein.
            if dist < options.interval:
                # Are they in the same frame? (and on the same strand)
                # If they are we assume it's an in-frame insertion
                if first_comp[2] == second_comp[2]:
                    TYPE = "INSERTION"
                # Are they on the same strand but a different frame?
                elif first_comp[2] * second_comp[2] > 0:
                    TYPE = "FRAMESHIFT"
                else:
                    TYPE = "INVERSION"
            if TYPE is not None:
                print "%s\t%s\t%d\t%d\t%d\t%d\t%s" %(contig, query, first_comp[0], first_comp[1], second_comp[0], second_comp[1], TYPE)

        # Check for stop codons (this probably doesn't strictly need to be done separately but...)
        for ii in range(len(ls)):
            TYPE = None
            contig_region = getContigRegion(whole_contig_seq, ls[ii][0], ls[ii][1], False, True)
            # Should I only be reporting these if they are not too close to the end?
            if "*" in contig_region:
                # Note - find will locate the first instance of "*" which is good... if there are multiple stop codons in a hit we're only interested in the first one.
                loc = contig_region.find("*") + 1
                if 100 - float(loc)/float(len(contig_region))*100.0 > options.nonsense_pct:
                    TYPE = "NONSENSE"
            if TYPE is not None:
                print "%s\t%s\t%d\t%d\t%d\t%d\t%s" %(contig, query, ls[ii][0], ls[ii][1], ls[ii][0], ls[ii][1], TYPE)
