#!/usr/bin/python


'''Test suite for tree functions'''
import unittest
from TreeFuncs import *
from ete2 import Tree

class RerootTreeTester(unittest.TestCase):
    def setUp(self):
        self.test_tree = Tree("test_tree.nwk")
        self.initial_nwk = self.test_tree.write(format=0)
        self.good_gene = "fig|2210.3.peg.1758"
        self.good_org = "2210.3"
    def test_bad_argments(self):
        self.assertRaises(ValueError, rerootEteTree, self.test_tree, root_leaf = "provided both leaf", root_leaf_part = "and a part of it")
        self.assertRaises(ValueError, rerootEteTree, self.test_tree, root_leaf = "this_is_not_in_the_tree")
        self.assertRaises(ValueError, rerootEteTree, self.test_tree, root_leaf_part = "neither_is_this")
    def test_good_arguments(self):
        # This gene is in the tree
        test_tree_rerooted = rerootEteTree(self.test_tree, root_leaf = self.good_gene)
        first_nwk = test_tree_rerooted.write(format=0)
        # So is this organism...
        test_tree_rerooted_2 = rerootEteTree(self.test_tree, root_leaf_part = self.good_org)
        second_nwk = test_tree_rerooted.write(format=0)
        # And they should be the same.
        self.assertEqual(first_nwk, second_nwk)
        # But should be different from the input tree
        self.assertNotEqual(self.initial_nwk, first_nwk)

class PretifyTreeTester(unittest.TestCase):
    def setUp(self):
        self.test_tree = Tree("test_tree.nwk")
        self.initial_nwk = self.test_tree.write(format=0)
        self.good_gene = "fig|2210.3.peg.1758"
        self.good_org = "2210.3"
    def test_prettify_tree(self):
        new_tree, new_ts = prettifyTree(self.test_tree)
    def test_standardize_tree(self):
        new_tree = standardizeTreeOrdering(self.test_tree)

if __name__ == "__main__":
    unittest.main()
