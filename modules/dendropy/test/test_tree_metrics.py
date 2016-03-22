#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests of tree metrics.
"""

import random
import unittest
import math
from cStringIO import StringIO

import dendropy
from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import extendedtest

class TreeUnaryMetricsTest(unittest.TestCase):

    def testNBar(self):
        trees = datagen.reference_tree_list()
        # trees = dendropy.TreeList.get_from_path(
        #         src=pathmap.tree_source_path("pythonidae.beast.mcmc.trees"),
        #         schema='nexus')
        expected_values = [
            7.818181818181818,
            7.515151515151516,
            7.666666666666667,
            8.727272727272727,
            8.757575757575758,
            8.636363636363637,
            8.727272727272727,
            8.757575757575758,
            8.727272727272727,
            8.727272727272727,
            8.575757575757576,
            ]
        for idx, tree in enumerate(trees):
            observed = tree.N_bar()
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_colless_tree_imbalance(self):
        trees = datagen.reference_tree_list()
        # for tree in trees:
        #     print tree.colless_tree_imbalance()
        expected_values = [
            0.3024193548387097,
            0.2540322580645161,
            0.2762096774193548,
            0.3548387096774194,
            0.35685483870967744,
            0.344758064516129,
            0.3548387096774194,
            0.35685483870967744,
            0.3548387096774194,
            0.3548387096774194,
            0.3407258064516129,
            ]
        for idx, tree in enumerate(trees):
            observed = tree.colless_tree_imbalance(normalize="max")
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_colless_tree_imbalance2(self):
        # library(apTreeshape)
        # data(hivtree.treeshape)
        # print(paste("colless, raw: ", colless(hivtree.treeshape), sep=""))
        # print(paste("colless, pda: ", colless(hivtree.treeshape, "pda"), sep=""))
        # print(paste("colless, yule: ", colless(hivtree.treeshape, "yule"), sep=""))
        # # [1] "colless, raw: 992"
        # # [1] "colless, pda: 0.369977836654251"
        # # [1] "colless, yule: 0.993137704712054"
        tree = dendropy.Tree.get_from_path(
                src=pathmap.tree_source_path("hiv1.nexus"),
                schema='nexus')
        self.assertAlmostEqual(tree.colless_tree_imbalance(normalize=None), 992)
        self.assertAlmostEqual(tree.colless_tree_imbalance(normalize="pda"), 0.3699778366542512686443)
        self.assertAlmostEqual(tree.colless_tree_imbalance(normalize="yule"), 0.9931377047120540924041)

    def test_sackin_index(self):
        # library(apTreeshape)
        # data(hivtree.treeshape)
        # print(paste("sackin, raw: ", sackin(hivtree.treeshape), sep=""))
        # print(paste("sackin, pda: ", sackin(hivtree.treeshape, "pda"), sep=""))
        # print(paste("sackin, yule: ", sackin(hivtree.treeshape, "yule"), sep=""))
        # # [1] "sackin, raw: 2028"
        # # [1] "sackin, pda: 0.756365980579457"
        # # [1] "sackin, yule: 0.822783440343329"
        tree = dendropy.Tree.get_from_path(
                src=pathmap.tree_source_path("hiv1.nexus"),
                schema='nexus')
        self.assertAlmostEqual(tree.sackin_index(normalize=None), 2028)
        self.assertAlmostEqual(tree.sackin_index(normalize="pda"), 0.756365980579457)
        self.assertAlmostEqual(tree.sackin_index(normalize="yule"), 0.822783440343329)

    def test_b1(self):
        trees = datagen.reference_tree_list()
        # for tree in trees:
        #     print tree.B1()
        expected_values = [
            18.686544011544008,
            16.862301587301587,
            18.012301587301586,
            15.803210678210679,
            15.803210678210679,
            16.219877344877347,
            15.80321067821068,
            15.80321067821068,
            15.803210678210679,
            15.80321067821068,
            16.10321067821068,
            ]
        for idx, tree in enumerate(trees):
            observed = tree.B1()
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_treeness(self):
        trees = datagen.reference_tree_list()
        # for tree in trees:
        #     print tree.treeness()
        expected_values = [
            0.82043976304486,
            0.30678033634423607,
            0.2686940663128338,
            0.2674702980152253,
            0.2731856127080352,
            0.26942308963183575,
            0.2764640737121644,
            0.26096444220828763,
            0.2846852453916621,
            0.2791363657987356,
            0.28304948441090816,
            ]
        for idx, tree in enumerate(trees):
            observed = tree.treeness()
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_gamma2(self):
        trees = datagen.reference_tree_list()
        # for tree in trees:
        #     print tree.pybus_harvey_gamma()
        expected_values = [
            6.690070011342222,
            -2.1016546214332665,
            -2.2071830302961493,
            -0.9868763184862083,
            -1.1223514055125514,
            -1.0914035287339103,
            -0.9432772103480326,
            -0.9855794349340775,
            -0.7566110136514949,
            -0.4693672063234924,
            0.08314644690264045,
            ]
        for idx, tree in enumerate(trees):
            observed = tree.pybus_harvey_gamma()
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def testPHGamma(self):
        newick_str = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        tree = dendropy.Tree(stream=StringIO(newick_str), schema="newick")
        self.assertAlmostEqual(tree.pybus_harvey_gamma(tree), 0.546276, 4)


class TreeCompareTests(extendedtest.ExtendedTestCase):

    def setUp(self):
        self.tree_list1 = datagen.reference_tree_list()
        self.tree_list2 = datagen.reference_tree_list()

    def testNonMutatingDistinctTaxonSetSameStructComparisons(self):
        tl1_ts = self.tree_list1.taxon_set
        tl2_ts = self.tree_list2.taxon_set
        self.assertIsNot(tl1_ts, tl2_ts)
        for i, t1 in enumerate(self.tree_list1):
            t2 = self.tree_list2[i]
            t1_ts = t1.taxon_set
            t2_ts = t2.taxon_set
            self.assertIsNot(t1_ts, t2_ts)
            self.assertEqual(t1.symmetric_difference(t2), 0)
            self.assertAlmostEqual(t1.euclidean_distance(t2), 0)
            self.assertAlmostEqual(t1.robinson_foulds_distance(t2), 0)
            self.assertIs(t1.taxon_set, t1_ts)
            self.assertIs(t2.taxon_set, t2_ts)
        self.assertIs(self.tree_list1.taxon_set, tl1_ts)
        self.assertIs(self.tree_list2.taxon_set, tl2_ts)

    def testSymmetricDifferences(self):
        expected = {
            (0,1):60, (0,2):60, (0,3):60, (0,4):60, (0,5):60, (0,6):60, (0,7):60, (0,8):60,
            (0,9):60, (0,10):60, (1,2):14, (1,3):24, (1,4):22, (1,5):20, (1,6):24, (1,7):22, (1,8):24,
            (1,9):22, (1,10):22, (2,3):18, (2,4):16, (2,5):16, (2,6):18, (2,7):16, (2,8):18, (2,9):16, (2,10):16,
            (3,4):4, (3,5):4, (3,6):0, (3,7):4, (3,8):0, (3,9):2, (3,10):4, (4,5):2, (4,6):4, (4,7):0, (4,8):4,
            (4,9):2, (4,10):4, (5,6):4, (5,7):2, (5,8):4, (5,9):2, (5,10):4, (6,7):4, (6,8):0, (6,9):2, (6,10):4,
            (7,8):4, (7,9):2, (7,10):4, (8,9):2, (8,10):4, (9,10):2,
        }
        for i, t1 in enumerate(self.tree_list1[:-1]):
            for j, t2 in enumerate(self.tree_list2[i+1:]):
                v = t1.symmetric_difference(t2)
                self.assertEqual(expected[(i, i+j+1)], v)
#                print "(%d,%d):%d," % (i, i+j+1, v),
#                if (i * i+j+1) % 6 == 0:
#                    print

    def testEuclideanDistances(self):
        expected = {
            (0,1):442.518379997, (0,2):458.269219125, (0,3):492.707662859, (0,4):457.731995932, (0,5):463.419798784, (0,6):462.181969494,
            (0,7):439.865064545, (0,8):462.3054297, (0,9):479.06569226, (0,10):544.720324057, (1,2):105.534825723, (1,3):168.86739068, (1,4):119.287056085, (1,5):127.221894919, (1,6):125.918517173,
            (1,7):102.290062347, (1,8):130.5296198, (1,9):154.336066685, (1,10):247.555999428, (2,3):89.1098950842, (2,4):45.5124918081,
            (2,5):52.2607244547, (2,6):53.0477320261, (2,7):62.1391636266, (2,8):59.898883066, (2,9):79.3921379438, (2,10):172.187021923,
            (3,4):73.4046806483, (3,5):61.7211889655, (3,6):63.308525227,
            (3,7):113.043429355, (3,8):64.9098905352, (3,9):43.9926843558, (3,10):91.395044825, (4,5):22.881252195, (4,6):24.686671743,
            (4,7):47.14854215, (4,8):30.4425119229, (4,9):58.4893274048, (4,10):158.948156946, (5,6):24.7029660833, (5,7):56.9022982438, (5,8):25.0745838358, (5,9):45.9638357231, (5,10):146.364107049,
            (6,7):56.1301333366, (6,8):20.3469798051, (6,9):43.429825221, (6,10):145.712937469, (7,8):58.1647873304, (7,9):89.4537113125, (7,10):197.098347126, (8,9):40.5187846693, (8,10):145.393476072,
            (9,10):111.210401924,
        }
        for i, t1 in enumerate(self.tree_list1[:-1]):
            for j, t2 in enumerate(self.tree_list2[i+1:]):
                v = t1.euclidean_distance(t2)
                self.assertAlmostEqual(expected[(i, i+j+1)], v)
#                print "(%d,%d):%s," % (i, i+j+1, v),
#                if (i * i+j+1) % 6 == 0:
#                    print

    def testRobinsonFouldsDistances(self):
        expected = {
            (0,1):1849.2928245, (0,2):2058.49072588, (0,3):2196.0995614, (0,4):1953.16064964, (0,5):1984.76411566, (0,6):1943.24487014,
            (0,7):1723.09194669, (0,8):1920.18504491, (0,9):1998.04696628, (0,10):2406.42091465, (1,2):508.212960297, (1,3):702.092000773, (1,4):579.45550447, (1,5):577.047914047, (1,6):596.881857714,
            (1,7):535.123132276, (1,8):611.28408319, (1,9):632.852687475, (1,10):857.759045631, (2,3):364.804588356, (2,4):283.907134148,
            (2,5):305.534136399, (2,6):318.128572842, (2,7):424.71989186, (2,8):351.751319705, (2,9):358.03680072, (2,10):531.731219604,
            (3,4):315.556017395, (3,5):271.016494089, (3,6):314.906668504,
            (3,7):517.444273417, (3,8):343.014958112, (3,9):266.498405531, (3,10):278.870282525, (4,5):133.642994808, (4,6):134.649689854,
            (4,7):260.010627711, (4,8):148.901405649, (4,9):187.954728978, (4,10):477.315325085, (5,6):150.970483084, (5,7):277.342311245, (5,8):144.704886539, (5,9):159.326519241, (5,10):449.629738145,
            (6,7):256.511541887, (6,8):119.487158128, (6,9):182.878241583, (6,10):493.201642403, (7,8):237.16728985, (7,9):296.353239488, (7,10):709.696300851, (8,9):171.021015022, (8,10):522.572965967,
            (9,10):435.439226227,
        }
        for i, t1 in enumerate(self.tree_list1[:-1]):
            for j, t2 in enumerate(self.tree_list2[i+1:]):
                v = t1.robinson_foulds_distance(t2)
                self.assertAlmostEqual(expected[(i, i+j+1)], v)
#                print "(%d,%d):%s," % (i, i+j+1, v),
#                if (i * i+j+1) % 6 == 0:
#                    print

class FrequencyOfSplitsTest(unittest.TestCase):

    def setUp(self):
        self.trees = dendropy.TreeList.get_from_path(
                src=pathmap.tree_source_path('pythonidae.random.bd0301.tre'),
                schema='nexus')

    def testCount1(self):
        split_leaves = ['Python regius', 'Apodora papuana']
        f = self.trees.frequency_of_split(labels=split_leaves)
        self.assertAlmostEqual(f, 0.04)

    def testRaisesIndexError(self):
        split_leaves = ['Bad Taxon', 'Apodora papuana']
        self.assertRaises(IndexError, self.trees.frequency_of_split, labels=split_leaves)

if __name__ == "__main__":
    unittest.main()
