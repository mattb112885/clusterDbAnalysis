#!/usr/bin/python

'''This file contains scripts for the purpose of locating key files and
directories in the installed package.

Note that these functions will ONLY work if this file is kept in the lib/ directory
relative to the root.'''

# I include these lines so that the lib/ folders are always in the
# PYTHONPATH
import sys
import os
sys.path.append(os.path.abspath(__file__))

def locateRootDirectory():
    fileloc = os.path.abspath(__file__)
    path = os.path.dirname(fileloc)
    rootdir = os.path.abspath(os.path.join(path, ".."))
#    rootdir = "/home/benedic1/ClusterDB_BACKUP_4_26_2013/"
    return rootdir

def locateDatabase(raiseError=True):
    rootdir = locateRootDirectory()
    dbpath = os.path.join(rootdir, "db", "DATABASE.sqlite")
    if not os.path.exists(dbpath) and raiseError:
        raise ValueError("ERROR: Database file not found in expected location %s\n" %(dbpath) )
    return dbpath

def locateOrganismFile(raiseError=True):
    rootdir = locateRootDirectory()
    orgpath = os.path.join(rootdir, "organisms")
    if not os.path.exists(orgpath) and raiseError:
        raise ValueError("ERROR: organisms file not found in expected location %s\n" %(orgpath) )
    return orgpath

def locateGroupsFile(raiseError=True):
    rootdir = locateRootDirectory()
    groupspath = os.path.join(rootdir, "groups")
    if not os.path.exists(groupspath) and raiseError:
        raise ValueError("ERROR: groups file not found in expected location %s" %(groupsspath) )
    return groupspath

