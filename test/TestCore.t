#!/usr/bin/python


'''Test suite for core cluster identification functions'''
import unittest
from CoreGeneFunctions import *

class CoreTester(unittest.TestCase):
    def setUp(self):
        self.orglist = [ "Methanosarcina acetivorans C2A" ]
        self.runid = "all_I_1.7_c_0.4_m_maxbit"
    def test_bad_argments(self):
        # Must have at least one option
        self.assertRaises(ValueError, findGenesByOrganismList, self.orglist, self.runid)
        # Contradictory options
        self.assertRaises(ValueError, findGenesByOrganismList, self.orglist, self.runid, all_org = True, none_org = True)

    def test_good_arguments(self):
        # Should not fail
        # I checked that this gives the same results as the existing script.
        # and also verified that they look at least to be correct.
        clusterlist = findGenesByOrganismList(self.orglist, self.runid, all_org = True, uniq_org = True, only_org = True)
        for c in clusterlist:
            print "%s\t%s" %(c[0], c[1])

if __name__ == "__main__":
    unittest.main()
