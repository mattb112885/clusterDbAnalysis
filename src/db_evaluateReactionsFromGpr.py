#!/usr/bin/python

import fileinput, optparse, sqlite3, sys, re
from locateDatabase import *

usage = "%prog -g [GPR_file] -i [RunID] > Rxn_presence_absence"
description = """This function takes a two-column table containing
gene protein reaction relationships (GPR) and produces a table
telling whether the reaction is present in each organism in a given cluster run
based on the presence\absence of individual genes from the clustering results.

The GPR_file should have exactly two columns:
- A column of reaction IDs (first column)
- A column of Gene-protein relationships [I.E. "GeneX and GeneY"] (second column)

This will only work if the gene IDs are the same in the GPRs as they are in the database.

The gene IDs MUST be formatted the same way as they are in the database, i.e. in the 
fig\|\d+\.\d+\.peg\.\d+ format
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gprfile", help="GPR file (required, D=None)", action="store", type="str", dest="gprfile", default=None)
parser.add_option("-i", "--runid", help="Run ID to use to identify presence\absence of genes (requried, D=None)", action="store", type="str", dest="runid", default=None)
parser.add_option("-o", "--or", help="Replace all AND in the input GPR with OR (useful for diagnosing issues with missing subunits). D = False, evaluate as written",
                  action="store_true", dest="repor", default=False)
(options, args) = parser.parse_args()

if options.gprfile is None:
    sys.stderr.write("ERROR: GPR file (-g) is a required argument to evaluateReactionsFromGpr.py \n")
    exit(2)

if options.runid is None:
    sys.stderr.write("ERROR: Run ID (-i) is a required argument to evaluateReactionsFromGpr.py \n")
    exit(2)

geneFinder = re.compile("fig\|\d+\.\d+\.peg\.\d+")

# Lets get a dictionary from reactions to GPRs first
rxn2gpr = {}
for line in open(options.gprfile, "r"):
    spl = line.strip("\r\n").split("\t")
    gpr = spl[1]
    if options.repor:
        gpr = gpr.replace("and", "or")
    rxn2gpr[spl[0]] = gpr

# Now lets get a list of genes
genelist = []
for rxn in rxn2gpr:
    genes = geneFinder.findall(rxn2gpr[rxn])
    if genes is None:
        sys.stderr.write("WARNING: No genes with expected format found in GPR %s" %(rxn2gpr[rxn]))
        continue
    for gene in genes:
        genelist.append(gene)

if len(genelist) == 0:
    sys.stderr.write("ERROR: No genes were found with expected formatting - dont forget to replace your IDs with those present in the database...\n")
    exit(2)

# Use the database to get lists of organisms in each cluster that contains a gene in the query GPR...
con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# We want to make sure we get all the organisms for a particular cluster... but we don't need
# all the genes for those.
query1 = "SELECT clusterid,geneid FROM clusterorgs WHERE clusterorgs.runid=? AND clusterorgs.geneid=?"
query2 = "SELECT organism FROM clusterorgs WHERE clusterorgs.runid=? AND clusterorgs.clusterid=?"
query3 = "SELECT DISTINCT organism FROM clusterorgs WHERE clusterorgs.runid = ?"
cluster2orgs = {}
cluster2genes = {}

orglist = set()
cur.execute(query3, (options.runid, ))
for res in cur:
    orglist.add(res[0])

for gene in genelist:
    cur.execute(query1, (options.runid, gene))
    # For genes we only care about the ones actually appearing in our GPRs.
    clusterid = None
    for res in cur:
        clusterid = str(res[0])
        geneid = str(res[1])
        if clusterid in cluster2genes:
            cluster2genes[clusterid].add(geneid)
        else:
            cluster2genes[clusterid] = set()
            cluster2genes[clusterid].add(geneid)
    # Now lets get what organisms are in that cluster.
    cur.execute(query2, (options.runid, clusterid))
    for res in cur:
        org = res[0]
        orglist.add(org)
        if clusterid in cluster2orgs:
            cluster2orgs[clusterid].add(org)
        else:
            cluster2orgs[clusterid] = set()
            cluster2orgs[clusterid].add(org)

con.close()

# For each organism we iterate over each of those clusters and see if it is in there.
# If not, we assign the genes in those clusters to FALSE and if they are to TRUE.
# Then we evaluate...
rxn2presence = {}
rxn2presence["orgs"] = orglist

syntaxerrors = set()
nameerrors = set()
badrxns = []
for org in orglist:
    gene2presenceabsence = {}
    for cluster in cluster2orgs:
        val = False
        if org in cluster2orgs[cluster]:
            val = True
        for gene in cluster2genes[cluster]:
            # To use EVAL we need to have variables (gene names) with valid python syntax.
            gen = gene.replace("|", "_").replace(".", "_")
            gene2presenceabsence[gen] = val
    # Now we SHOULD be able to do eval.
    # If we have issues it means there's a formatting error in the GPR - we should print the offending ones and continue on.
    for rxn in rxn2gpr:
        gpr = rxn2gpr[rxn]
        # Replace in the gpr as well.
        badgpr = False
        try:
            rxnPresent = eval(gpr.replace("|", "_").replace(".", "_"), {"__builtins__":None}, gene2presenceabsence)
        except SyntaxError:
            syntaxerrors.add("WARNING: There was a syntax error (probably a missing parenthesis) in the following GPR: \n%s\t%s \n" %(rxn, gpr))
            badrxns.append(rxn)
            break
        except NameError:
            nameerrors.add("WARNING: The following GPR had a NameError (likely caused by a bad gene name in the GPR): \n%s\t%s \n " %(rxn,gpr))
            badrxns.append(rxn)
            break
        if rxn in rxn2presence:
            rxn2presence[rxn].append(rxnPresent)
        else:
            rxn2presence[rxn] = [rxnPresent]

for s in syntaxerrors:
    sys.stderr.write(s)
for s in nameerrors:
    sys.stderr.write(s)

for rxn in rxn2presence:
    # I do it this kind of round-about way because eval isn't always consistent with whether or not it throws a syntax error!
    # This makes me a bit worried about the accuracy of the results...
    if rxn in badrxns:
        continue
    print "%s\t%s" %(rxn, "\t".join( [ str(a) for a in rxn2presence[rxn] ] ))
