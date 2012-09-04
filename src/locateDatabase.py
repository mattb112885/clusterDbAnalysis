#!/usr/bin/python

import os

# This function will work as long as this file is in the same source directory as the rest of the source files...
def locateDatabase():
    fileloc = os.path.abspath(__file__)
    path = os.path.dirname(fileloc)
    s = os.path.abspath(os.path.join(path, "..", "db", "DATABASE.sqlite"))
    return s

