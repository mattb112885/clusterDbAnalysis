#!/usr/bin/python

# Calcualte score based on defined methods
#
# -c: Cutoff
# -m: Method (must match one of the methods below)
#
# -m minbit: Bit score / min( query self-bit, target self-bit)
# -m maxbit: Bit score / max( query self-bit, target self-bit)
# -m avgbit: Bit score * 2 / (query self-bit + target self-bit)

import optparse

validmethods = ['minbit', 'maxbit', 'avgbit']
