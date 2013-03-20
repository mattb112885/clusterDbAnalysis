#!/usr/bin/python

# This is a pipe command. Pipe in the results from the blastres_bitscore table
# that you want to calculate scores for.
#
# (e.g. from the specific organisms blast getter)
#
# Calcualte score based on defined methods
#
# -c: Cutoff
# -m: Method (must match one of the methods below)
#
# -m minbit: Bit score / min( query self-bit, target self-bit)
# -m maxbit: Bit score / max( query self-bit, target self-bit)
# -m maxselfbit: Bit score / query self-bit
# -m avgbit: Bit score * 2 / (query self-bit + target self-bit)
# 
# Local alignments (if you WANT domain hits)
# -m normhsp: Bit score / hsp length (equivalent to r in mcxdeblast)

import optparse, fileinput, sys

from FileLocator import *
from ClusterFuncs import *

symmetric_methods, unsymmetric_methods = getValidBlastScoreMethods()

usage="%prog -m [method] -c [cutoff] [options] < blast_result_table"
description = "Given a blast score table (augmented with self-bit scores for query and target genes), calculates a similarity value based on the desired scoring metric. Currently implemented metrics: %s" %(" ".join(symmetric_methods))
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--method", help="Method name", action="store", type="str", dest="method", default=None)
parser.add_option("-c", "--cutoff", help="Score cutoff to use", action="store", type="float", dest="cutoff", default=None)
parser.add_option("-n", "--noprint", help="Set this flag if you do not want to print edges less than the cutoff as zeros (default = False - print those edges)", action="store_true", dest="noprint", default=False)
(options, args) = parser.parse_args()

if options.method not in symmetric_methods:
    sys.stderr.write("ERROR: Method %s not in list of valid (symmetric) methods: %s\n" %(" ".join(symmetric_methods)) )
    exit(2)

# Note - I limit to 100000 at a time due to possible memory issues dealing with
# huge amounts of BLAST results to compute.
blastres = []
maxn = 100000
n=0
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    blastres.append(spl)
    n+=1
    if n >= maxn:
        scorelist = calculateScoreFromBlastres(blastres, options.method, options.cutoff, include_zeros = not options.noprint, needsymmetric = True)
        for score in scorelist:
            print "%s\t%s\t%s" %(score[0], score[1], str(score[2]))
        n = 0
        blastres = []

scorelist = calculateScoreFromBlastres(blastres, options.method, options.cutoff, include_zeros=not options.noprint, needsymmetric = True)
for score in scorelist:
    print "%s\t%s\t%s" %(score[0], score[1], str(score[2]))
