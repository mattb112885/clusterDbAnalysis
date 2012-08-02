#!/usr/bin/python

# EMACS does not support tabs very well - trying to insert a tab often ends up inserting
# a bunch of spaces instead.
#
# We therefore need a way to make a tab-delimted row to insert into e.g. the groups
# or organisms files
#
# I could use VIM but instead I wrote this silly script.
#
# Input arguments: Any number of things you want separated by tabs
#
# Output arguments: Tab-delimited those things
#
# Make sure anything with special characters in it is surrounded in single quotes when calling from bash
# so it is interpreted literally...

import optparse, fileinput

usage="""%prog [Things to separate with tabs] > Tab-delimited line
%prog -i [Input file] > Tab_delimited line"""
description = "Turn all arguments into a single row separated by tabs"

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-i", "--infile", help="Combine all lines of the input file and separate them by tabs. Use \"-\" for input from stdin", action="store", type="str", dest="infile", default=None)
(options, args) = parser.parse_args()

st = ""
if options.infile == None:
    st = "\t".join(args[0:])
elif options.infile == "-":
    st = "\t".join( [ s.strip("\r\n") for s in fileinput.input("-") ] )
else:
    st = "\t".join( [ s.strip("\r\n") for s in open(options.infile, "r") ] )

print st
