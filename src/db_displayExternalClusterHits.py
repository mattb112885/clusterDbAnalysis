#!/usr/bin/env python

# FIXME - Put these into a library (should be their own library since external clusters
# are kind of unique (definitely distinct from de novo computed ones)
def getHitsToExternalClusters(external_clusterid_list, evalue, cur):
    '''
    This function identifies all genes in the database that have a RPSBLAST hit with an e-value better
    than the specified evalue (closer to 0) and returns a list tuples (geneID, external_clusterID, E-value)

    You must run setup_step4.sh for this function to work.

    external_clusterid_list is a list of strings (e.g. pfam00001) corresponding to external cluster IDs.
    evalue is a floating-point value (e.g. 1E-5)
    cur is a SQLite cursor object to the ITEP database.

    The function returns a unique list of gene IDs (if there are hits to multiple input clusters the gene ID
    is only returned once).
    '''

    hits = []
    q = """SELECT querygene,evalue FROM rpsblast_results
         INNER JOIN external_clusters ON external_clusters.cdd_id = rpsblast_results.cdd_id
         WHERE external_clusters.external_clusterid = ? AND evalue < ?"""
    for externalid in external_clusterid_list:
        cur.execute(q, (externalid, evalue))
        for res in cur:
            hits.append( (res[0], externalid, res[1]) )
    return hits

def getRpsBlastForQueryGenes(querygene_list, evalue, cur, database = "all"):
    '''
    Given a list of query gene IDs (querygene_list), returns all available RPSBLAST hits to external databases with
    specified significance level.

    database: If database is selected to be other than "all" (e.g. PFAM), the function only returns results that hit that
    database rather than everything.
    cur: A SQLite cursor object pointing to the ITEP database

    Returns a list of tuples corresponding to RPSBLAST hit table rows.
    '''
    
    addl = ""
    if options.database != "all":
        lk = "\"%s%%\"" %(options.database)
        addl = " AND external_clusters.external_clusterid LIKE %s" %(lk)

    cmd = """SELECT rpsblast_results.*, external_clusters.external_clusterid FROM rpsblast_results
         INNER JOIN external_clusters ON external_clusters.cdd_id = rpsblast_results.cdd_id
         WHERE rpsblast_results.querygene = ? AND evalue < ? %s """ %(addl)

    hits = []
    for geneid in querygene_list:
        cur.execute(cmd, (geneid, evalue) )
        for res in cur:
            hits.append( tuple(res) )
    return hits

import fileinput
import operator
import optparse
import os
import sqlite3
import sys

from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors as rcolor

from BioPythonGraphics import *
from ClusterFuncs import *
from FileLocator import *

usage = "%prog [options] < gene_id_list"
description = """ Use Biopython to display the hits to external clusters
(via RPSBLAST) mapped onto a protein or a list of proteins. You have to run setup_step4.sh
before this function will work.
"""

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-d", "--database", help="Only display hits for this database (e.g. PFAM, COG, cd, ...) (D: Show all of them in separate rows above the gene)",
                  action="store", type="str", dest="database", default="all")
parser.add_option("-e", "--evalue", help="E-value cutoff for hits to display (D:1E-5)",
                  action="store", type="float", dest="evalue", default=1E-5)
parser.add_option("-g", "--geneocl", help="Column number for gene ID, starting from 1 (D=1)",
                  action = "store", type="int", dest="gc", default=1)
parser.add_option("-o", "--outdir", help="Output directory for all of the PNG files for the input genes (D: externalHitGraphics)",
                  action="store", type="str", dest="outdir", default="externalHitGraphics")
parser.add_option("-s", "--showevalue", help="Show E-value for RPSBLAST hits along with the names (D: Show names only)",
                  action="store_true", dest="showevalue", default=False)
parser.add_option("-m", "--maxhits", help="Maximum number of hits to display on the figure (D:Not limited). Try this if the figure gets too messy.",
                  action="store", type="int", dest="maxhits", default=None)
(options, args) = parser.parse_args()

if options.outdir is None:
    sys.stderr.write("ERROR: output directory (-o) is a required argument\n")
    exit(2)

# Make it not fail if the directory doesn't exist already.
if not os.path.isdir(options.outdir):
    if os.path.exists(options.outdir):
        sys.stderr.write("ERROR: Specified output folder name already exists as an ordinary file\n")
        exit(2)
    os.mkdir(options.outdir)

# Read input
gc = options.gc - 1
genelist = set()
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    genelist.add(spl[gc])

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Set up a BioPython display object for the gene and its RPSBLAST hits, and display them.
for geneid in genelist:
    # Get RPSBLAST hit results for this gene
    rpsblast_hits = getRpsBlastForQueryGenes( [geneid], options.evalue, cur, database=options.database )

    # Get a divergent color mapping.
    rps_targets = list(map(operator.itemgetter(12), rpsblast_hits))
    colorscheme = colormap(rps_targets)

    # RPSBLAST hits have coordinates relative to the beginning of the protein (not relative to the location on the contig).
    # This is good, since it makes our indexing much less of a headache.
    # But it also means we need to get the length of the gene.
    geneinfo = getGeneInfo( [geneid], cur )
    start = 1
    end = len(geneinfo[0][11]) + 1
    annotation = geneinfo[0][9]

    gd_diagram = GenomeDiagram.Diagram(geneid)
    gd_track_for_gene = gd_diagram.new_track(1, name=geneid)
    gd_feature_set = gd_track_for_gene.new_set()

    # Some basic properties of the figure itself
    arrowshaft_height = 0.3
    arrowhead_length = 0.3
    scale = 100     #AA per px for the diagram
    bordercol=rcolor.white

    color = rcolor.black

    gene_feature = SeqFeature(FeatureLocation(start, end), strand=+1, id = geneid)

    gd_feature_set.add_feature(gene_feature, name = gene_feature.id + "(%s)" %(annotation),
                               color=color, border = bordercol,
                               sigil="ARROW", arrowshaft_height=arrowshaft_height, arrowhead_length = arrowhead_length,
                               label=True, label_position = "start", label_angle = 0, label_size = 24
                               )

    n = 2
    for rps in rpsblast_hits:
        rps_id = rps[1]
        hitstart = rps[6]
        hitend = rps[7]
        evalue = rps[10]
        rps_name = rps[12]
        color = colorscheme[rps_name]
        
        # FIXME - this should become its own function
        q = "SELECT clustername FROM external_clusters WHERE cdd_id = ?"
        cur.execute(q, (rps_id,))
        for res in cur:
            cname = str(res[0])

        rps_name = rps_name + "(%s)" %(cname)
        if options.showevalue:
            rps_name += " - Evalue=%s" %(evalue)

        gd_track = gd_diagram.new_track(n, name=rps_name)
        gd_feature_set = gd_track.new_set()
        hit_feature = SeqFeature( FeatureLocation(hitstart, hitend), strand=+1, id=rps_name)
        gd_feature_set.add_feature(hit_feature, name = rps_name,
                                   color=color, border = bordercol,
                                   label=True, label_position = "start", label_angle = 0, label_size = 24
                                   )
        n += 1
        if options.maxhits is not None and n > (options.maxhits + 1):
            break

    gd_diagram.draw(format="linear", start=start, end=end, fragments=1)
    pathname = os.path.join( options.outdir, sanitizeString(geneid, False) + ".svg" )
    gd_diagram.write(pathname, "SVG")
    sys.stderr.write("Image file for gene %s written to %s\n" %(geneid, pathname))

con.close()
