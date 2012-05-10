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
# -m avgbit: Bit score * 2 / (query self-bit + target self-bit)
# 
# Local alignments (if you WANT domain hits)
# -m normhsp: Bit score / hsp length (equivalent to r in mcxdeblast)

def minbit(bitscore, qselfbit, tselfbit):
    return float(bitscore)/min( float(qselfbit), float(tselfbit) )

def maxbit(bitscore, qselfbit, tselfbit):
    return float(bitscore)/max( float(qselfbit), float(tselfbit) )

def avgbit(bitscore, qselfbit, tselfbit):
    return float(bitscore) * 2 / (float(qselfbit) + float(tselfbit) )

def normhsp(bitscore, hsplen):
    return float(bitscore)/float(hsplen)

import optparse, fileinput

validmethods = ['minbit', 'maxbit', 'avgbit', 'normhsp']

parser = optparse.OptionParser()
parser.add_option("-m", "--method", help="Method name", action="store", type="str", dest="method", default=None)
parser.add_option("-c", "--cutoff", help="Score cutoff to use", action="store", type="float", dest="cutoff", default=None)
(options, args) = parser.parse_args()

if not options.method in validmethods:
    print "ERROR: Invalid method passed to db_makeBlastScoreTable."
    print "Passed method:"
    print options.method
    print "Current valid methods are:"
    print validmethods
    exit(2)

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
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
    elif options.method == 'avgbit':
        score = avgbit(bitscore, qselfbit, tselfbit)
    elif options.method == 'normhsp':
        score = normhsp(bitscore, hsplen)

    if score < options.cutoff:
        score = 0

    print "%s\t%s\t%1.3f" % (qgene, tgene, score)
