#!/usr/bin/python

# Split space-delimited string into separate lines
# useful for piping...
import fileinput

for line in fileinput.input("-"):
    spl = line.strip().split(" ")
    for s in spl:
        print s
    
