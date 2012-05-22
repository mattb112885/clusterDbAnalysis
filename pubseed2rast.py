#!/usr/bin/python

# This is a pipe command intended to be called as a subscript to
# convertToRast.sh
#
# Convert a PUBSEED table with annotations attached to it
# into RAST format so it can be uploaded into the database.
#

import fileinput

for line in fileinput.input("-"):
    spl = line.strip().split("\t")

    # Coding sequences only.
    if not spl[6] == "CDS":
        continue

    # In the order of what we want from the output...
    contig = spl[2]
    featureid = spl[0]
    featuretype = "peg"
    # We don't use this for anything so lets just keep it as the contig ID
    location = contig
    genestart = spl[3]
    geneend = spl[4]
    genestrand = spl[5]
    annotation = spl[1]
    aliases = ""
    figfam = ""
    evidence = spl[11]
    ntseq = spl[13] # These must be added by calling svr_fasta (like the sh script does) before calling this function
    aaseq = spl[14]

    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, featureid, featuretype, location, genestart, geneend, genestrand, annotation, aliases, figfam, evidence, ntseq, aaseq)
    
