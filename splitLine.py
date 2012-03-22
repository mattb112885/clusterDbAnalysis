#!/usr/bin/python

import fileinput

for line in fileinput.input("-"):
    spl = line.strip().split(" ")
    for s in spl:
        print s
    
