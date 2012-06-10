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

import sys

print "\t".join(sys.argv[1:])
