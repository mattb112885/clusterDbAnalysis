#!/usr/bin/python


'''Test suite for core cluster identification functions'''
import unittest
from CoreGeneFunctions import *
from TreeFuncs import *
from ete2 import Tree

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

class AddCoreTester(unittest.TestCase):
    def setUp(self):
        self.runid = "methanosarcinales_Methanocella_I_1.7_c_0.3_m_maxbit"
        self.ete_tree = Tree("test_org_tree.nwk")
        self.ete_tree = rerootEteTree(self.ete_tree, root_leaf = "Methanocella_paludicola_SANAE")
        self.ete_tree, self.ts = prettifyTree(self.ete_tree, show_bootstraps = False)
    def testCore(self):
        (newtree, clusterdata) = addCoreDataToTree(self.ete_tree, self.runid, all_org = True, color="Black")
        # Add second one (that was the whole point of having it be a function)
        (newtree, clusterdata) = addCoreDataToTree(newtree, self.runid, all_org = True, only_org = True, color="Red")
        # Make sure both of them show up.
        newtree.show(tree_style = self.ts)

if __name__ == "__main__":
    unittest.main()
