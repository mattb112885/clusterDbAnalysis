#!/usr/bin/python

import optparse, os, re, sys
from Bio import SeqIO
from locateDatabase import *

usage = """%prog -l mysql_loginname -p mysql_password -d mysql_database_string [options]
%prog -f orthomcl_config_file [options] """
description = """
  WARNING - This script is not yet complete, do NOT try to do analysis based on it yet!!

  This is a wrapper script for converting our data into a format in which it can be run with OrthoMCL.
  I separated this from the wrapper for actually RUNNING OrthoMCL so that we can run multiple times with different settings
  but the same BLAST results.

  The script requires installation of MYSQL and having a database created.  It also requires having the orthomcl binaries
  (in $ORTHOMCLROOT/bin) added to your PATH variable. See OrthoMCL help for details on the steps performed here."""

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-l", "--login", help="Login name for MySQL", action="store", type="str", dest="login", default=None)
parser.add_option("-p", "--password", help="Password for MySQL (required)", action="store", type="str", dest="password", default=None)
parser.add_option("-d", "--dbstring", 
                  help=""" Database string for the MySQL database (required). 
                           For a central install it will look like this:
                           dbi:MySql:[database_name]
                           For a local install it will instead look like:
                           dbi:MySql:[database_name]:localhost:[port]""",
                  action="store", type="str", dest="dbstring", default=None)
parser.add_option("-o", "--orthotable", help="Name of ortholog table to create in database (D: Ortholog)", action="store", type="str", dest="orthotable", default="Ortholog")
parser.add_option("-i", "--inparalogtable", help="Name of inparalog table to create in database (D: InParalog)", action="store", type="str", dest="inparalogtable", default="InParalog")
parser.add_option("-t", "--coortholog", help="Name of co-ortholog table to create in database (D: CoOrtholog)", action="store", type="str", dest="coorthologtable", default="CoOrtholog")
parser.add_option("-e", "--evalcut", help="Cutoff for log-E value (int, D=-5, i.e. E-values higher than E-5 are treated as 0)", action="store", type="int", dest="evalcut", default=-5)
parser.add_option("-c", "--pctcut", help="Percent match cutoff for OrthoMCL (int, D=50)", action="store", type="int", dest="pctcut", default=50)
parser.add_option("-n", "--newconfigfile", help="Name of new config file to create (D: orthomcl.new.config)", action="store", type="str", dest="newconfigfile", default="orthomcl.new.config")
parser.add_option("-f", "--configfile", help="A previously-generated orthoMCL config file (D: Provide -l, -p, and -d instead and we will create a file)",
                  action="store", type="str", dest="configfile", default=None)

(options, args) = parser.parse_args()

# Either -l, -p, and -d must all be specified OR a previously-generated config file must be provided.
if options.configfile is None and (options.login is None or options.password is None or options.dbstring is None):
    sys.stderr.write("""ERROR: Either a previously-existing OrthoMCL config file must be specified (-f), 
                        or a login, password, and database connection string must be specified by command line (-l, -p, and -d) 
                        so one can be generated\n""")
    exit(2)

# The alternative here is to have the command-line options override the config file? I'm only comfortable doing this
# if the config file has an empty field there to begin with.
if options.configfile is not None and (options.login is not None or options.password is not None or options.dbstring is not None):
    sys.stderr.write("""WARNING: Multiple sets of logins, passwords, and/or dbstrings provided - 
                        will use the existing config file and ignore command-line arguments!\n""")

#####################################
# STEPS FOR ORTHOMCL Running        #
# 1) Install MySQL [assumed done]   #
#####################################

#####################################
# 2) Install MCL [assumed done]     #
#####################################

#####################################
# 3) Configure OrthoMCL             #
#####################################

# If there isn't a config file passed in already, create one.
if options.configfile is None:
    fid = open(options.newconfigfile, "w")
    fid.write("dbVendor=mysql\n")
    fid.write("dbConnectString=%s\n" %(options.dbstring))
    fid.write("dbLogin=%s\n" %(options.login))
    fid.write("dbPassword=%s\n" %(options.password))
    # For these I just use the defaults
    fid.write("similarSequencesTable=SimilarSequences\n")
    fid.write("orthologTable=Ortholog\n")
    fid.write("inParalogTable=InParalog\n")
    fid.write("coOrthologTable=CoOrtholog\n")
    # I think there is only one possible value for this one.
    fid.write("interTaxonMatchView=InterTaxonMatch\n")
    # Cutoffs
    fid.write("percentMatchCutoff=%d\n" %(options.pctcut))
    fid.write("evalueExponentCutoff=%d\n" %(options.evalcut))
    fid.write("oracleIndexTblSpc=NONE\n")
    fid.close()
    sys.stderr.write("New config file written to %s\n" %(options.newconfigfile))
    options.configfile = options.newconfigfile

#####################################
# 4) Run orthomclInstallSchema      #
#####################################

# Have orthoMCL run SQL to set up the database schema...
# We don't necessarily want to re-load all the BLAST info (takes a long time) if
# it hasn't changed at all.
sys.stderr.write("Installing OrthoMCL schema...\n")
os.system("orthomclInstallSchema %s" %(options.configfile))

########################################
# 5-7) - Make BLAST results files in   #
# the required format.                 #
# The required format is as follows:   #
# [organismid]|[geneid]                #
# (the gene ID should be sanitized to  #
# not have pipes in it)                #
#                                      #
# Also need to make FASTA file with    #
# that format.                         #
########################################

# The BLAST parser requires all of the organisms to be in their own fasta file
# and the fasta files must be named "organismid".fasta
# Also... apparently the regex that must be matched is
# $fastaFile =~ /(\w+).fasta/
# \w+ means alphanumeric plus "_"
# This means we have to sanitize the periods out of the fasta file names and the organism IDs.
#
# Also note that though the "sql" file indicates that the maximum size of gene IDs is 15 characters the
# perl script that actually is used to make the schema uses a maximum size of 40. 40 is enough for us.

sys.stderr.write("Generating a FASTA file with orthomcl-compliant formatting...\n")

# I make the fasta in a separate folder because orthomcl asks for a directory name
# and searches that for all the fasta files (and yells at you if anything in there doesn't
# have the right format).
orthofastadir = os.path.join(os.path.dirname(locateDatabase()), "orthofasta")

origfastadir = os.path.join(os.path.dirname(locateDatabase()), "..", "faa")
faafinder = re.compile(".*\.faa$")
organismIdFinder = re.compile("\d+\.\d+")

for f in os.listdir(origfastadir):
    # If it is a FASTA (FAA) file
    if faafinder.search(f) is not None:
        origfasta = os.path.join(origfastadir, f)
        recs = SeqIO.parse(open(origfasta, "r"), "fasta")
        # Modify the sequence IDs to conform to orthoMCL format.
        ids2seqs = {}
        for rec in recs:
            myid = rec.id
            myseq = rec.seq
            myid = myid.replace("fig|", "")
            match = organismIdFinder.search(myid)
            if match is None:
                sys.stderr.write("ERROR: Gene ID %s not in expected format so I cannot get the organism ID out\n" %(myid))
                exit(2)
            orgid = match.group(0)
            # Sanitize out periods in the organism ID
            orgid = orgid.replace(".", "_")
            # OrthoMCL requires organism ID to be separated from the gene ID by a pipe.
            myid = "%s|%s" %(orgid, myid)
            ids2seqs[myid] = myseq
        outfile = open(os.path.join(orthofastadir, orgid + ".fasta"), "w")
        for myid in ids2seqs:
            outfile.write(">%s\n%s\n" %(myid, ids2seqs[myid]))
        outfile.close()

# ... and we need to ensure that the genes in the BLAST results have the same format.
# This will take some time but WAY less than re-running BLAST.
sys.stderr.write("Reformatting BLAST output file to use orthoMCL-compliant IDs...\n")
blastresfile = os.path.join(os.path.dirname(locateDatabase()), "blastres_cat")
orthoblastres = os.path.join(os.path.dirname(locateDatabase()), "blastres_cat_orthomcl")

try:
    outfile = open(orthoblastres, "r")
    outfile.close()
except IOError:
    outfile = open(orthoblastres, "w")
    for line in open(blastresfile, "r"):
        spl = line.strip("\r\n").split("\t")
        qid = spl[0]
        qid = qid.replace("fig|", "")
        qorg = organismIdFinder.search(qid)
        qid = "%s|%s" %(qorg.group(0).replace(".", "_"), qid)
        
        tid = spl[1]
        tid = tid.replace("fig|", "")
        torg = organismIdFinder.search(tid)
        tid = "%s|%s" %(torg.group(0).replace(".", "_"), tid)
        
        spl[0] = qid
        spl[1] = tid
        outfile.write("\t".join(spl) + "\n")
    outfile.close()

########################################
# 8) run orthomclBlastParser - it adds #
# necessary info abot query proteins   #
# to the table.                        #
########################################

# Im' HOPING this works.
sys.stderr.write("Running orthomcl Blast parser to put blast results in the correct format for output...\n")
orthomclBlastParserFile = os.path.join(os.path.dirname(locateDatabase()), "blastres_orthomcl_modfiied")
try:
    fid = open(orthomclBlastParserFile, "r")
    fid.close()
except IOError:
    cmd = "orthomclBlastParser %s %s > %s" %(orthoblastres, orthofastadir, orthomclBlastParserFile)
    os.system(cmd)

########################################
# 9) orthomclLoadBlast - load up the   #
#  mySQL database                      #
########################################

sys.stderr.write("Reloading the MySQL database with re-formatted BLAST output...(this will take a long time)\n")
cmd = "orthomclLoadBlast %s %s" %(options.configfile, orthomclBlastParserFile)
os.system(cmd)

sys.stderr.write("Done with loading. Run orthoMclWrapper.py to complete the rest of the steps with desired cluster settings\n")
