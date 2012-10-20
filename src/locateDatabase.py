#!/usr/bin/python

import os
import optparse
import sys

def main():
    description = """This file contains a function for internal use - locating the sqlite database to allow the scripts
to be run from anywhere without explicitly specifying the location in a config script"""
    parser = optparse.OptionParser(description=description)
    parser.parse_args()
    return 0

if __name__ == "__main__":
    main()

# This function will work as long as this file is in the same source directory as the rest of the source files...
def locateDatabase():
    fileloc = os.path.abspath(__file__)
    path = os.path.dirname(fileloc)
    s = os.path.abspath(os.path.join(path, "..", "db", "DATABASE.sqlite"))
    return s

