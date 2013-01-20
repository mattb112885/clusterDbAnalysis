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
    return rootdir

def locateDatabase():
    rootdir = locateRootDirectory()
    dbpath = os.path.join(rootdir, "db", "DATABASE.sqlite")
    if not os.path.exists(dbpath):
        sys.stderr.write("ERROR: Database file not found in expected location %s" %(dbpath) )
        raise ValueError
    return dbpath

def locateOrganismFile():
    rootdir = locateRootDirectory()
    orgpath = os.path.join(rootdir, "organisms")
    if not os.path.exists(orgpath):
        sys.stderr.write("ERROR: organisms file not found in expected location %s" %(orgpath) )
        raise ValueError
    return orgpath

def locateGroupsFile():
    rootdir = locateRootDirectory()
    groupspath = os.path.join(rootdir, "groups")
    if not os.path.exists(groupspath):
        sys.stderr.write("ERROR: groups file not found in expected location %s" %(groupsspath) )
        raise ValueError
    return groupspath

