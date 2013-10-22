#!/usr/bin/env python

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

def sanitizeByType(container, sanitizeby='tsv', onlycolumns=False):
    '''for a iterable of strings, carry out sanitizeString by:
        line, 
        tsv (all or onlycolumns), 
        fasta headers, or 
        leaf in nwk'''
    
    assert sanitizeby in set(['line', 'tsv', 'newick', 'fasta'])
    if sanitizeby=='line': 
        for line in container:
            print sanitizeString(line.strip("\r\n"), False)
    if sanitizeby=='tsv': 
        for line in container:
            if onlycolumns: 
                newline = line.strip("\r\n").split("\t")
                for i in onlycolumns: 
                    newline[i-1]=sanitizeString(newline[i-1], False)
            else:
                newline=[sanitizeString(item.strip("\r\n"), False) for item in line.split("\t")]
            print "\t".join(newline)
    if sanitizeby=='newick':
        from ete2 import Tree
        t=Tree("".join(container))
        for l in t:
            l.name=sanitizeString(l.name, False)
        print t.write()
    if sanitizeby=='fasta': 
        from Bio import SeqIO
        from StringIO import StringIO
        from sys import stdout
        fasta = StringIO("".join(container))
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seq_record.id=sanitizeString(seq_record.description, False)
            seq_record.description=''
            SeqIO.write(seq_record, stdout, "fasta")
