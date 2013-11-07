#!/usr/bin/env python

from FileLocator import *
import os, optparse, sys, glob

usage = "%prog [options] [\"searchstring_1\" \"searchstring_2\" ... ] > filelist"
description="""List all files in the ITEP directories provided as part of this software package.
If searchstrings are provided, only returns those matching the search strings."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-w", "--maxw", 
                  help="Maximum number of characters wide (D=Print each file on its own line)", 
                  action="store", type="int", dest="maxw", default=1)
#can be used for auto-documentation like: 
#db_listDbFiles.py -s  | xargs -I % bash -c "echo;echo "-------------------------------";echo %; % -h" > DOCS.txt
parser.add_option("-s", "--simpleoutput", 
                  help="Force printing exactly one command per line, no other formatting, overides -w (D=Print human readable output)", 
                  action="store_true", dest="simpleoutput", default=False)
(options, args) = parser.parse_args()

subpaths = ['src','scripts','src/utilities']
paths = [os.path.join(locateRootDirectory(), sub) for sub in subpaths]

def printcol(ls, indent=0, out=sys.stdout, width=72):
    maxw = width-indent
    # Pretty-print into columns, but keep alphabetical order by column
    maxl = max( [ len(s) for s in ls ] )+1 #need 1 space to prevent 2 longest items abutting
    gte_one = lambda x: x if x>1 else 1
    cols = gte_one(int(maxw//maxl)) #floor division gives us the number of columns, must be >1, leave space for indent
    rows = (len(ls) + (cols-1)) // cols #this is a formula for ceiling division
    padded = sorted(ls)+[""]*(len(ls)%cols) #add entries
    justified = [s.ljust(maxl) for s in padded] #add spaces to entries so they are the same length
    lines = [justified[r::rows] for r in range(rows)] # make alphabetized by column
    text = "\n".join([" "*indent+"".join(l) for l in lines])
    out.write(text)
    out.write("\n") #need final return

for path in paths:
    ls = glob.glob(os.path.join(path,"*")) # will not list directories or .hidden, gives path
    ls = [f for f in ls if (os.access(os.path.join(path,f),os.X_OK) and os.path.isfile(f))] #only executable files
    ls = [os.path.basename(f) for f in ls] #get filenames
    query= set([q.lower() for q in args])
    def hits(h): return any([h.lower().find(q)>-1 for q in query])
    if args==[]:
        found=ls
    else:
        found = [h for h in ls if hits(h)]
    #keep only if it's executable
    found=[h for h in found]
    if len(found) > 0:
        if options.simpleoutput: 
            printcol(found, width=1)
        else:
            print('\nPrograms found in '+path+':')
            printcol(found, indent=2, width = options.maxw)
    elif not options.simpleoutput: 
        print('\nNo programs found in '+path+'.')
