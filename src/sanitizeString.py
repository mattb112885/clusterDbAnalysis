#!/usr/bin/python

# Standard function for string sanitization
# Mostly to prevent Newick files from getting confused

import re, sys

def sanitizeString(string, warnOfReplacement):
    rem = re.compile("[^\.0-9A-Za-z]")
    s = rem.sub("_", string)
    if warnOfReplacement and not string == s:
        sys.stderr.write("WARNING: String %s replaced with sanitized version  %s\n" %(string, s))
    return s
