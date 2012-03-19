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

import optparse, fileinput

validmethods = ['minbit', 'maxbit', 'avgbit']

parser = optparse.OptionParser()
parser.add_option("-m", "--method", help="Method name", action="store", type="str", dest="method", default=None)
parser.add_option("-c", "--cutoff", help="Score cutoff to use", action="store", type="float", dest="method", default=None)
(options, args) = parser.parse_args()

if not options.method in validmethods:
    print "ERROR: Invalid method passed to db_makeBlastScoreTable."
    print "Current valid methods are:"
    print validmethods
    exit(2)

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    # Required information
    bitscore = spl[
