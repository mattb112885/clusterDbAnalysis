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

def minbit(bitscore, qselfbit, tselfbit):
    return float(bitscore)/min( float(qselfbit), float(tselfbit) )

def maxbit(bitscore, qselfbit, tselfbit):
    return float(bitscore)/max( float(qselfbit), float(tselfbit) )

def maxselfbit(bitscore, qselfbit):
    return float(bitscore)/float(qselfbit)

def avgbit(bitscore, qselfbit, tselfbit):
    return float(bitscore) * 2 / (float(qselfbit) + float(tselfbit) )

def normhsp(bitscore, hsplen):
    return float(bitscore)/float(hsplen)

import optparse, fileinput
from locateDatabase import *

validmethods = ['minbit', 'maxbit', 'maxselfbit', 'avgbit', 'normhsp']

usage="%prog -m [method] -c [cutoff] [options] < blast_result_table"
description = "Given a blast score table (augmented with self-bit scores for query and target genes), calculates a similarity value based on the desired scoring metric. Currently implemented metrics: %s" %(" ".join(validmethods))
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--method", help="Method name", action="store", type="str", dest="method", default=None)
parser.add_option("-c", "--cutoff", help="Score cutoff to use", action="store", type="float", dest="cutoff", default=None)
parser.add_option("-n", "--noprint", help="Set this flag if you do not want to print edges less than the cutoff as zeros (default = False - print those edges", action="store_true", dest="noprint", default=False)
(options, args) = parser.parse_args()

if not options.method in validmethods:
    sys.stderr.write("ERROR: Invalid method passed to db_makeBlastScoreTable\n")
    sys.stderr.write("Current valid methods are:\n%s\n" %("\t".join(validmethods)))
    exit(2)

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    # Required information
    bitscore = spl[11]
    qselfbit = spl[12]
    tselfbit = spl[13]
    hsplen = spl[3]

    qgene = spl[0]
    tgene = spl[1]

    if options.method == 'minbit':
        score = minbit(bitscore, qselfbit, tselfbit)
    elif options.method == 'maxbit':
        score = maxbit(bitscore, qselfbit, tselfbit)
    elif options.method == 'maxselfbit':
        score = maxselfbit(bitscore, qselfbit)
    elif options.method == 'avgbit':
        score = avgbit(bitscore, qselfbit, tselfbit)
    elif options.method == 'normhsp':
        score = normhsp(bitscore, hsplen)

    if score < options.cutoff:
        if options.noprint:
            continue
        else:
            score = 0

    print "%s\t%s\t%1.3f" % (qgene, tgene, score)
