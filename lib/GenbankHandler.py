#!/usr/bin/python

'''
This file contains functions intended to deal with modifications to GENBANK files
for increased compatibility between ITEP and other tools.

'''

import sys

def addItepGeneIdsToGenbank(multi_gbk_object, tbl, truncateContigIds=False):
    '''
    Given a RAW table (list of tuples) and a genbank object,
    attempt to match up every object in the RAW table with every
    gene in the table and add the ITEP IDs as a db_xrer (called ITEP)
    to the genbank file.

    If truncateContigIds is set to TRUE the contig IDs for every locus are truncated
    to 16 characters.
    '''
    # Index search on the table
    startidx = 4
    stopidx = 5
    seqidx = 11
    ididx = 1
    indexed_array = {}
    for ln in tbl:
        # Do both of these because of possible incompatibilities in how the orderinf of
        # start and stop works...
        start = int(ln[startidx])
        stop = int(ln[stopidx])
        seq = ln[seqidx]
        myid = ln[ididx]
        indexed_array[(start, stop)] = (seq, myid)
        indexed_array[(stop, start)] = (seq, myid)

    # First we need this to not be a generator so we can actually modify the darn thing
    multi_gbk_object = list(multi_gbk_object)

    # I match up by the following:
    # Start location
    # Stop location
    # DNA sequence (due to possible differences in translation)
    # If these all match with a given element of the table then its good.
    all_contigs = set()
    for ii in range(len(multi_gbk_object)):
        # If we have to, truncate the contig names so that biopython can actually produce some output.
        # However, if the user doesn't tell us to truncate them and they would have to, throw an error instead
        # to warn the user of the problem.
        if len(multi_gbk_object[ii].name) > 16:
            if truncateContigIds:
                sys.stderr.write("WARNING: Truncated contig ID %s to %s (16 characters)\n" %(multi_gbk_object[ii].name, multi_gbk_object[ii].name[:16]))
                multi_gbk_object[ii].name = multi_gbk_object[ii].name[:16]
            else:
                raise IOError("ERROR: Contig ID %s is too long - writing it back out with biopython will fail" %(multi_gbk_object[ii].name))
            pass

        # Truncating at 16 characters could cause the same contig ID to appear mutliple times in the same file, which is BAD.
        # Lets check for such nonsense...
        if multi_gbk_object[ii].name in all_contigs:
            raise IOError("ERROR: Contig name %s (before or after truncation) would have appeared multiple times in the same genbank file. This is bad!")
        else:
            all_contigs.add(multi_gbk_object[ii].name)
            pass

        # Add ITEP IDs to the genbank files.
        for jj in range(len(multi_gbk_object[ii].features)):
            # We don't want to modify things that aren't coding sequences...
            if multi_gbk_object[ii].features[jj].type != "CDS":
                continue
            record = multi_gbk_object[ii].features[jj].extract(multi_gbk_object[ii])
            seq = record.seq.tostring()
            # Get locations to match up (note they are splice indexes and start from 0 hence the +1 to start and nothing added to end)
            location = multi_gbk_object[ii].features[jj].location
            featurestart = int(location.start + 1)
            featureend = int(location.end)
            querytup = (featurestart, featureend)
            if querytup in indexed_array and seq == indexed_array[querytup][0]:
                if "db_xref" in multi_gbk_object[ii].features[jj].qualifiers:
                    multi_gbk_object[ii].features[jj].qualifiers["db_xref"].append("ITEP:%s" %(indexed_array[querytup][1]))
                else:
                    multi_gbk_object[ii].features[jj].qualifiers["db_xref"] = [ "ITEP:%s" %(indexed_array[querytup][1]) ]
            else:
                sys.stderr.write("BAD - the gene in the following location did not have a match in the table file!\n")
                sys.stderr.write(str(location) + "\n")
                sys.stderr.write(str(seq) + "\n")
            pass
        pass

    return multi_gbk_object
