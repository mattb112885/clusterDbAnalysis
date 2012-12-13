#! /usr/bin/env python
# -*- coding: utf-8 -*-


#This is a pipe command
#
#Given a set of taxonids (from stdin) it will download all genebank files for that taxonid into files named with the taxonid
#
#Includes various functions to get genbank files from NCBI and search between databases at NCBI
#
#More information avaliavble by running with the -h flag.  

#this script uses Entrez searching.  Note that there are interlinks between genome, taxonomy, and nuccore, but that we are able to search using one of these two forms, as everything can be limited by taxonomy.
# also note the serch syntax is different for different databases
# taxonomy can be filtered by: "taxonomy genome"[Filter] AND txid402880[Organism:noexp]
# nuccore can be filtered by: txid39152[Organism:noexp] AND srcdb_refseq[Properties]
# see http://www.ncbi.nlm.nih.gov/books/NBK21100/#A286 for some details about the terms and limits

from Bio import Entrez, SeqIO
import os, sys, re, getpass, socket
import argparse, textwrap 
from collections import defaultdict

def uniq(seq):
   '''Order-preserving unique function'''
   noDupes = []
   [noDupes.append(i) for i in seq if not noDupes.count(i)]
   return noDupes

def getGenbanks(idlist, outlocation='', overwrite=True):
    '''Download genbank files using a unique identifier (ID, GID, or nucleotide accession) '''
    gbfiles = []
    for gbid in idlist:
        filename = getGenbank(gbid, outlocation=outlocation, overwite=overwrite)
        #store filelocation in dict
        gbfiles.append(filename)
    return gbfiles

def getGenbank(gbid, outlocation='', overwrite=True):
    '''Download a single genbank file (auxiliary to getGenbanks)'''
    filename = os.path.join(outlocation, str(gbid) +".gbk")
    out_handle = open(filename, "w")
    fileexists = not os.path.isfile(filename)
    #get gb file from NCBI
    net_handle = Entrez.efetch(db="nucleotide",id=gbid, rettype="gb", retmode="text")
    if overwrite or fileexists :
        out_handle.write(net_handle.read())
        sys.stderr.write("\t Genbank for %s downloaded and saved at %s\n" % (gbid, filename))
    else: sys.stderr.write("WARNING %s file exists and/or overwrite not selected\n" % filename)
    out_handle.close()
    net_handle.close()
    return filename

def strainID_to_genbankIDs(strainID, refseqonly = True):
    '''Get all nucliotide sequences with this taxonID (chromosomes, plasmids, etc.). 
    Note that the taxon must be the lowest designator to get only one organism's data'''
    nucIDs=[]
    term = "txid%s[Organism:noexp] AND nuccore assembly[filter]" % strainID
    #we may only want refseqs - note that this is a different limit than other dbs (like genome = refseq[filter])
    if refseqonly == True:
        term = term + " AND srcdb_refseq[Properties]"
    genomeIDs = Entrez.read(Entrez.esearch(db="nucleotide", term = term))
    nucIDs = nucIDs + genomeIDs['IdList']
    #make unique (should be already)
    nucIDs = list(set(nucIDs))
    sys.stderr.write("Strain %s is nucleotide ID(s) %s\n" % (strainID, nucIDs))
    return nucIDs

def accession_to_strainID(accession, return_nucID = True):
    nucID  = accessionID_to_genbankIDs(accession)
    strainIDs=[]
    summaries = Entrez.read(Entrez.esummary(db="nucleotide", id=nucID))
    for s in summaries:
        strainIDs.append(s['TaxId'])
    strainIDs = list(set(strainIDs))
    sys.stderr.write("Accession %s is taxonID %s\n" % (accession, strainIDs))
    assert len(strainIDs)==1, "WARNING: more than one ID returned from query of NCBI %s" % accession
    if return_nucID==True:
        return strainIDs[0], nucID
    else:
        return strainIDs[0]

def accessionID_to_genbankIDs(accession):
    '''Search NCBI for the ID (the unique numerical key)'''
    genbank = Entrez.read(Entrez.esearch(db="nucleotide", term=accession))
    IDlist = genbank['IdList']
    #TODO: find the most recent version?
    sys.stderr.write("Accession %s is nucleotide ID  %s\n" % (accession, IDlist))
    assert len(IDlist) == 1, "WARNING: more than one NCBI ID returned for an accesntion number"
    return IDlist

def taxonIDs_to_strainIDs(taxonIDs, hasgenome=True, donotexpand=False, excludewgs=False):
    '''Find all individual strains in a taxonID'''
    strainIDs = []
    for taxonID in taxonIDs:
        #see http://www.ncbi.nlm.nih.gov/books/NBK21100/#A286 for some details about the terms and limits
        query = 'txid%s[subtree] AND "taxonomy genome"[filter] AND terminal[prop]' % taxonID
        if hasgenome:
            query = query + ' AND (("taxonomy assembly"[Filter]) OR ("taxonomy genome"[Filter]) OR ("taxonomy genome2"[Filter]))'
        if excludewgs:
            query = query + ' AND NOT wgs[filter]'
        NCBIdata = Entrez.read(Entrez.esearch(db="taxonomy", term=query))
        #there can me multiple records returned
        taxon_strainIDs = NCBIdata['IdList']
        if len(taxon_strainIDs)<1 :
            sys.stderr.write("WARNING: No sequence information of the type you requested for %s\n" % taxonID)
        if donotexpand:
            assert len(taxon_strainIDs)==1, "WARNING: more than one strain ID returned from query of NCBI %s" % taxonID
        else:
            if len(taxon_strainIDs)>1:
                #only need to let user know if there is more than one, otherwise the numbers are the same
                sys.stderr.write("Note: multiple strains associated with taxonID %s: %s\n" % (taxonID, taxon_strainIDs))
        strainIDs = strainIDs + taxon_strainIDs
    return strainIDs

def gb_taxon(gbfile, try_bioparsing=False):
    '''Get a NCBI taxonomy ID from a standard or RAST genbank file, looking in several places
    try_bioparsing will fail on many files, try at your own risk'''
    #open genbank or multi-genbank
    gb_handle = open(gbfile,"r")
    #scan for regex (from Matt's methods in /convertNcbiFilesToRast.pl), note that db_xref, Db_xref, dbxref, and others are all found in GB files.  Sigh.
    taxIdFinder = re.compile("xref=.*taxon:(\d*)")
    for line in gb_handle:
        # Discard comment lines
        if line.startswith("#"):
            continue
        try: strainID = taxIdFinder.search(line).group(1)
        except AttributeError: continue
        else:
            sys.stderr.write("File %s is taxonID %s\n" % (gbfile, strainID))
            return strainID
    gb_handle.close()
    if try_bioparsing==True:
        #try with parser
        gb_seqrec = SeqIO.parse(gb_handle, "genbank")
        for sequence in gb_seqrec:
            #get from entire locus dbxref if it is there
            if sequence.dbxrefs:
                try: dbxrefs = dict([xref.split(':') for xref in sequence.dbxrefs])
                except ValueError: pass
                else:
                    try: return dbxrefs['taxon']
                    except KeyError: pass
            #get from feature if it is not
            for feat in sequence.features: # source should be first one, but this will loop
                try: dbxrefs = dict([xref.split(':') for xref in feat.qualifiers['db_xref']])
                except ValueError: pass
                else:
                    try: return dbxrefs['taxon']
                    except KeyError: pass
    #if a taxonID is not in the file, perhapse it has an accession?   so, look it up in NCBI
    accession = gb_accession(gbfile)
    strainID = accession_to_strainID(accession, return_nucID = False)
    sys.stderr.write("File %s has accession %s identified as taxonID %s by NCBI\n" % (gbfile, accession, strainID))
    return strainID

def gb_accession(gbfile):
    '''Get a NCBI accession ID from a standard or RAST genbank file'''
    #open genbank or multi-genbank
    gb_handle = open(gbfile,"r")
    gb_seqrec = SeqIO.parse(gb_handle, "genbank")
    sequence = gb_seqrec.next()
    accession = sequence.name
    gb_handle.close()
    sys.stderr.write("File %s has accession %s\n" % (gbfile, accession))
    return accession

def NCBIversion(db):
    record = Entrez.read(Entrez.einfo(db=db), validate=False)
    update = record["DbInfo"]["LastUpdate"]
    sys.stderr.write("This script is using the current NCBI %s database, version %s\n" % (db, update))

def NCBIdbinfo(db):
    '''Helper function for developers to look up NCBI DB structure'''
    record = Entrez.read(Entrez.einfo(db=db), validate=False)
    record["DbInfo"]['Description']
    info = record["DbInfo"]
    printinfo = [(x['Name'], x['Description']) for x in info['FieldList']]
    for item in printinfo:
        print "%s: %s \n" % item

def concatinate_files(dictionary, outlocation = '', extension="gbk"):
    '''Given a dictionary from an output filename to a list of input filenames,
    concatinates the input files into the output file.
    The provided extension is added to the output file name.

    file1 --> [infile1, infile2]
    file2 --> [infile2, infile3]

    would create two output files, file1.extension (concatinating infile1 and infile2)
    and file2.extension (concatinating infile2 and infile3).
    '''
    for filenamebase, files in dictionary.items():
        outfilename = os.path.join(outlocation, str(filenamebase) + "." + extension)
        outfile = open(outfilename, 'w')
        for f in files:
            outfile.write(open(f).read())
            sys.stderr.write("\t concatinating %s into %s\n" % (f, outfilename))
        outfile.close()


#Entrez requires an email, if not set namually, this guesses one
email = None
guessemail = getpass.getuser() + '@' + socket.getfqdn()
if email == None:
    email = guessemail
Entrez.email = email

if __name__=="__main__":

    usage='''
    examples:
        Download_genbank.py -t -s 406327 418008 523844
        OR,  with a file of 1 taxonID per line named taxonlist.txt:
        cat taxonlist.txt | Download_genbank.py -t
        OR
        Download_genbank.py -t taxonlist.txt
    '''
    description='''
    Given a set of taxonids, NCBI nucleotide accessions, or genbank files (from stdin) this script will 
download all genebank files into one file named with the taxonid. Note that for files and 
accessions, only one genbank file will be downloaded. For taxonids, all sub-taxons will 
be downloaded.
    '''
    #arguments and help
    __version__ = '0.0.2'
    parser = argparse.ArgumentParser(description=description, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent(usage)
                                     )

    parser.add_argument('-v', '--version', action='version', version=__version__)
    inputstrains = parser.add_mutually_exclusive_group(required=True)
    inputstrains.add_argument("-a", "--accessions",
                              help="input is a list of NCBI accession numbers",
                              action="store_true", dest="accessioninput")
    inputstrains.add_argument("-t", "--taxonIDs",
                              help="input is a list of taxonID numbers",
                              action="store_true", dest="taxonidinput")
    inputstrains.add_argument("-g", "--genbank-files",
                              help="input is a list of genbank files",
                              action="store_true", dest="genbankinput")
    parser.add_argument("-s", "--string-input",
                              help="override input, use this string",
                              action="store", metavar = "LIST_OF_IDS",
                              nargs='+', dest="stringinput")
    parser.add_argument("-d", "--outdir",
                        help = "Output directory for download (D: /tmp)",
                        action = "store", dest="outputdir",
                        default = "/tmp")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin, help="without a file, uses stdin")
    inargs = parser.parse_args()

    outlocation = ''
    #report where we are getting data
    NCBIversion('nucleotide')
    NCBIversion('taxonomy')
    sys.stderr.write("\n")

    #first, where is the input coming from?
    #tests: NC_017527 NC_007796 NC_007955 or 523846, 868131, 402880 or '/tmp/56708791.gbk', '246200.gbk'
    try: len(inargs.stringinput) > 0
    except TypeError:
      inputlist = inargs.infile
      #if it came from a file, will need to chomp the newlines
      inputlist=[line.strip() for line in inputlist]
    else:
      inputlist = inargs.stringinput

    #convert from input list to what we want to download
    #strainIDlist are a list of straions, for each one there will be a list of genbank nucleotide IDs
    if inargs.accessioninput:
        #these are all 1:1
        IDs = [accession_to_strainID(accession, return_nucID = True) for accession in inputlist]
        strainIDlist, genbankIDs = zip(*IDs)
    if inargs.genbankinput:
        #these are all 1:1 with the input, but we need to look up the genbankID from the accestions in a seperate step
        IDs = [(gb_taxon(gbfile), gb_accession(gbfile)) for gbfile in inputlist]
        strainIDlist, accessionIDs = zip(*IDs)
        genbankIDs = [accessionID_to_genbankIDs(accession) for accession in accessionIDs]
    if inargs.taxonidinput:
        #these are not 1:1 with the input if we are recursing (or if there are multiple refseqs per organism
        strainIDlist = uniq(taxonIDs_to_strainIDs(inputlist, hasgenome = True))
        genbankIDs = [strainID_to_genbankIDs(strainID, refseqonly = True) for strainID in strainIDlist]

    #make a dict for each organism ID, it will hold lists of genbank ids to download, that will be replaced with filenames
    #use a list defualt so we can append all the nucliotide IDs that have the same organism ID so they end up in the same file
    organismIDs = defaultdict(list)
    for st, gb in zip(strainIDlist, genbankIDs):
        organismIDs[st] = organismIDs[st] + gb
    #get one or more gb file for each entry (organism OR actual genbank ID)
    for strainID, genbankIDs in organismIDs.items():
        if len(genbankIDs ) > 0:
            sys.stderr.write("Data for strain %s\n" % strainID)
            gbfiles = [getGenbank(genbankID, outlocation=inargs.outputdir, overwrite=True) for genbankID in genbankIDs]
            organismIDs[strainID] = gbfiles
            #concatinate into a file, only getting one accession if that was the input.
            concatinate_files(organismIDs)
            # Clean up.
            for f in gbfiles:
               sys.stderr.write("Removing intermediate file %s\n" %(f))
               os.remove(f)

    #tests:
    '''
    #
    #simple 2 accession case
    runfile(r'/home/jamesrh/Documents/code/Download_genbank.py', args=r'-a -s NC_007955 NC_017527', wdir=r'/home/jamesrh/Documents/code')
    #simple file case
    runfile(r'/home/jamesrh/Documents/code/Download_genbank.py', args=r'-g -s /tmp/91772082.gbk /tmp/386000717.gbk', wdir=r'/home/jamesrh/Documents/code')
    #simple taxonID case (this should pick up a chromosome and a plasmid)
    runfile(r'/home/jamesrh/Documents/code/Download_genbank.py', args=r'-t -s 259564 1110509', wdir=r'/home/jamesrh/Documents/code')
    #recursive taxonID case (this should pick up chromosome and plasmid for at least one organism)
    runfile(r'/home/jamesrh/Documents/code/Download_genbank.py', args=r'-t -s -o 143067', wdir=r'/home/jamesrh/Documents/code')
    '''

