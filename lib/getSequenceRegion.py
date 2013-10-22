#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

######
# getContigRegion(seq, start, stop, wantProtein)
#
# For DNA sequences,
# use start < stop for + strand and start > stop for - strand
#
# Use 1-based counting (start=1 means the first nucleotide in the contig
#   and end = 3 means stop on the third nucleotide)
#
# - strand sequences are automatically reverse-complemented.
#
# If wantProtein is specified, the function translates the sequence
# and returns that instead.

def getContigRegion(myseq, start, stop, isprot, translate):
    if stop < start and isprot:
        raise IOError("If a protein sequence is specified, the start location must be less than the stop")
    if isprot and translate:
        raise IOError("Cannot ask for translation for a protein sequence")

    if start < stop:
        initialSequence = myseq[start-1:stop]
    else:
        initialSequence = myseq[stop-1:start]

    seqobj = Seq(initialSequence)
    if stop < start and not isprot:
        seqobj = seqobj.reverse_complement()

    if translate:
        seqobj = seqobj.translate()

    # Return the final sequence.
    return str(seqobj)
