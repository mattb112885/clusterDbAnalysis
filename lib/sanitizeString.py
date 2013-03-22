#!/usr/bin/python

# Standard function for string sanitization
# Mostly to prevent Newick files from getting confused

import re, warnings

'''This is an internal function to standardize sanitation of strings for unput into various formats.
Any character that is not a letter or number is replaced with an underscore.'''

def sanitizeString(string, warnOfReplacement):
    # I would've kept periods in here but RangerDTL chokes on them...
    # sigh. That's all I can say.
    rem = re.compile("[^0-9A-Za-z]")
    s = rem.sub("_", string)
    if warnOfReplacement and not string == s:
        warnings.warn("WARNING: String %s replaced with sanitized version  %s\n" %(string, s) )
    return s

