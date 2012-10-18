#!/usr/bin/python

import optparse, os, re, sys
from Bio import SeqIO
from locateDatabase import *

usage = """%prog -l mysql_loginname -p mysql_password -d mysql_database_string [options]
%prog -f orthomcl_config_file [options] """
description = """
  WARNING - This script is not yet complete, do NOT try to do analysis based on it yet!!

  This is a wrapper script for converting our data into a format in which it can be run with OrthoMCL and then
  running it with the specified settings.

  The script requires installation of MYSQL and having a database created.  It also requires having the orthomcl binaries
  (in $ORTHOMCLROOT/bin) added to your PATH variable. See OrthoMCL help for details on the steps performed here.

  If redundant settings are present in the config file and in the inputs, the settings in the config file are overridden by
  the command line input. """

parser = optparse.OptionParser(usage=usage, description=description)

# Configuration options
parser.add_option("-f", "--configfile", help="A previously-generated orthoMCL config file (D: Provide -l, -p, and -d instead and we will create a file)",
                  action="store", type="str", dest="configfile", default=None)
parser.add_option("-n", "--newconfigfile", help="Name of new config file to create with command-line settings (D: orthomcl.new.config)", 
                  action="store", type="str", dest="newconfigfile", default="orthomcl.new.config")
parser.add_option("-l", "--login", help="Login name for MySQL", action="store", type="str", dest="login", default=None)
parser.add_option("-p", "--password", help="Password for MySQL (required)", action="store", type="str", dest="password", default=None)
parser.add_option("-d", "--dbstring", 
                  help=""" Database string for the MySQL database (required). 
                           For a central install it will look like this:
                           dbi:MySql:[database_name]
                           For a local install it will instead look like:
                           dbi:MySql:[database_name]:localhost:[port]""",
                  action="store", type="str", dest="dbstring", default=None)
# Overwritable options for the config file
parser.add_option("-o", "--orthotable", help="Name of ortholog table to create in database (D: Ortholog)", action="store", type="str", dest="orthotable", default="Ortholog")
parser.add_option("-i", "--inparalogtable", help="Name of inparalog table to create in database (D: InParalog)", action="store", type="str", dest="inparalogtable", default="InParalog")
parser.add_option("-t", "--coortholog", help="Name of co-ortholog table to create in database (D: CoOrtholog)", action="store", type="str", dest="coorthologtable", default="CoOrtholog")
parser.add_option("-e", "--evalcut", help="Cutoff for log-E value (int, D=-5 i.e. cut off if E-value is bigger than 1E-5)", action="store", type="int", dest="evalcut", default=-5)
parser.add_option("-c", "--pctcut", help="Percent match cutoff for OrthoMCL (int, D=50, i.e. cutoff if percent match is less than 50)", action="store", type="int", dest="pctcut", default=None)
# Other options for other MCL programs
parser.add_option("-a", "--inflation", help="Inflation value for MCL (D: 1.5)", action="store", type="float", dest="inflation", default=1.5)
parser.add_option("-g", "--logfile", help="Orthomcl pairs Log file (D: Make a dummy name and delete it afterwards)", action="store", type="str", dest="logfile", default=None)
parser.add_option("-r", "--forcereload", help="Force reload of the database with BLAST results (D: Only re-load if the reformatted BLAST info file is newly created)",
                  action="store_true", dest="forcereload", default=False)
parser.add_option("-k", "--keeptemp", help="Keep temporary files made by orthoMCL (D: Delete them - except the new config file which is stored in the filename specified by -n)",
                  action="store_true", dest="keeptemp", default=False)

(options, args) = parser.parse_args()

# Either -l, -p, and -d must all be specified OR a previously-generated config file must be provided.
if options.configfile is None and (options.login is None or options.password is None or options.dbstring is None):
    sys.stderr.write("""ERROR: Either a previously-existing OrthoMCL config file must be specified (-f), 
                        or at least a login, password, and database connection string must be specified by command line (-l, -p, and -d) 
                        so one can be generated\n""")
    exit(2)

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

def getWantedOption(line, optionfinder, cmdlineoption):
    match = optionfinder.search(line)
    option = match.group(1)
    value = match.group(2)
    if len(value) == 0:
        if cmdlineoption is None:
            sys.stderr.write("""ERROR: No option for configuration parameter %s
            found on command line or in the orthomcl config file.\n""" %(option))
            raise IOError
        else:
            return cmdlineoption
    else:
        if cmdlineoption is None:
            return value
        elif value == str(cmdlineoption):
            return value
        else:
            sys.stderr.write("WARNING: Command line value for parameter %s will overwrite value in orthomcl config file\n" %(option))
            return cmdlineoption

# If there is a config file passed in, take the values from that but overwrite them with command-line arguments.
if options.configfile is not None:
    optionfinder = re.compile("^(.*)=(.*)$")
    for line in open(options.configfile, "r"):
        if line.startswith("#"):
            continue
        line = line.strip("\r\n")
        match = optionfinder.search(line)
        if match is None:
            continue

        if match.group(1) == "dbVendor":
            continue
        elif match.group(1) == "dbConnectString":
            options.dbstring = getWantedOption(line, optionfinder, options.dbstring)
        elif match.group(1) == "dbLogin":
            options.login = getWantedOption(line, optionfinder, options.login)
        elif match.group(1) == "dbPassword":
            options.password = getWantedOption(line, optionfinder, options.password)
        elif match.group(1) == "similarSequencesTable":
            continue
        elif match.group(1) == "orthologTable":
            options.orthotable = getWantedOption(line, optionfinder, options.orthotable)
        elif match.group(1) == "inParalogTable":
            options.inparalogtable = getWantedOption(line, optionfinder, options.inparalogtable)
        elif match.group(1) == "coOrthologTable":
            options.coorthologtable = getWantedOption(line, optionfinder, options.coorthologtable)
        elif match.group(1) == "interTaxonMatchView":
            continue
        elif match.group(1) == "percentMatchCutoff":
            options.pctcut = int(getWantedOption(line, optionfinder, options.pctcut))
        elif match.group(1) == "evalueExponentCutoff":
            options.evalcut = int(getWantedOption(line, optionfinder, options.evalcut))
        elif match.group(1) == "oracleIndexTblSpc":
            continue
        else:
            sys.stderr.write("WARNING: Option %s in config file not recognized so it will be skipped\n" %(match.group(1)))
            continue
        
# Create a new config file containing all the user-specified info and
# whatever defaults were there for what is not user-specified
fid = open(options.newconfigfile, "w")
fid.write("dbVendor=mysql\n")
fid.write("dbConnectString=%s\n" %(options.dbstring))
fid.write("dbLogin=%s\n" %(options.login))
fid.write("dbPassword=%s\n" %(options.password))
# For these I just use the defaults
fid.write("similarSequencesTable=SimilarSequences\n")
fid.write("orthologTable=%s\n" %(options.orthotable))
fid.write("inParalogTable=%s\n" %(options.inparalogtable))
fid.write("coOrthologTable=%s\n" %(options.coorthologtable))
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
    sys.stderr.write("Reformatted BLAST file %s already exists\n" %(orthoblastres))
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

sys.stderr.write("Running orthomcl Blast parser to put blast results in the correct format for output...\n")
orthomclBlastParserFile = os.path.join(os.path.dirname(locateDatabase()), "blastres_orthomcl_modfiied")
try:
    fid = open(orthomclBlastParserFile, "r")
    fid.close()
    sys.stderr.write("orthomcl Blast parser results file already exists\n")
    reloadDb = False
except IOError:
    cmd = "orthomclBlastParser %s %s > %s" %(orthoblastres, orthofastadir, orthomclBlastParserFile)
    os.system(cmd)
    reloadDb = True

########################################
# 9) orthomclLoadBlast - load up the   #
#  mySQL database                      #
########################################

if reloadDb or options.forcereload:
    sys.stderr.write("Reloading the MySQL database with re-formatted BLAST output...(this will take a long time)\n")
    cmd = "orthomclLoadBlast %s %s" %(options.configfile, orthomclBlastParserFile)
    os.system(cmd)
else:
    sys.stderr.write("BLAST results file in correct format already existed before running this script so it was assumed that the DB was loaded already. If this is not the case specify --forcereload\n")
    
if options.logfile is None:
    options.logfile = "ORTHOMCL_LOGFILE"

# Check if the output file for the specified settings already exists.
# Otherwise why bother re-doing all of this?

expectedfilename = "orthomcl_all_I_%1.4f_c_%d_%d" %(options.inflation, options.evalcut, options.pctcut)
expectedoutput = os.path.join(os.path.dirname(locateDatabase()), "..", "clusters", expectedfilename)
if not (options.forcereload or reloadDb):
    try:
        open(expectedoutput, "r")
        sys.stderr.write("Expected output file %s already exists...\n" %(expectedoutput))
        exit(0)
    except IOError:
        sys.stderr.write("OrthoMCL output %s does not already exist so creating it..\n" %(expectedoutput))

# Now that sanity checks are over lets run orthoMCL. First get pairs...
# cleanup="yes" is necessary to make it possible to run orthoMCL multiple times with different settings.
# Otherwise, it will refuse to run because a table that it attempts to create already exists.
cmd1 = "orthomclPairs %s %s cleanup=%s" %(options.configfile, options.logfile, "yes")
sys.stderr.write("Calculating orthoMCL pairs files... (this will take a long time)\n")
os.system(cmd1)

# Now we have pairs and we can extract them...
cmd2 = "orthomclDumpPairsFiles %s" %(options.configfile)
sys.stderr.write("Dumping orthoMCL results to files...\n")
os.system(cmd2)

# ... and then run MCL on the resulting "MclInput" file
# I make a temporary output file because I will need to re-make the gene IDs
# in the format expected by our database.
tmp_mcl_output = "mclOut_BADIDS"

cmd3 = "mcl mclInput --abc -I %1.4f -o %s" %(options.inflation, tmp_mcl_output)
sys.stderr.write("Running MCL on the orthoMCL outputs... (this could take some time)\n")
os.system(cmd3)

# Now I make the output file in the expected format.
# This regex MUST be lazy because there are multiple genes on a line
orgReplacer = re.compile("\S+?\|")
fout = open(expectedoutput, "w")
for line in open(tmp_mcl_output, "r"):
    line = orgReplacer.sub("fig|", line)
    fout.write(line)
fout.close()

# Finally, remove the temporary files created by orthoMCL
# NOT INCLUDING the new config file!!!
if not options.keeptemp:
    os.remove(options.logfile)
    os.remove(tmp_mcl_output)
    os.remove("mclInput")
    # rmdir only works on empty directories.
    os.system("rm -r pairs")
    
