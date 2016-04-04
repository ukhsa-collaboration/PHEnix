'''
Created on 8 Dec 2015

@author: alex
'''
import unittest

import vcf

import os.path as path
from phe.variant import VariantSet
from phe.variant_filters import PHEFilterBase
from phe.variant_filters.DepthFilter import DepthFilter


class Test(unittest.TestCase):

    def setUp(self):
        filter_config = {"min_depth":5, "ad_ratio":0.9}
        self.filter = DepthFilter(filter_config)

        vcf_file = path.join(path.dirname(__file__), "sample.vcf")

        self.var_set = VariantSet(vcf_file, filter_config)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def test_consensus(self):
        expected_consensus = ["W"] + 8 * ["N"] + ['G', 'W', 'Y', 'N', 'S', 'R', 'C']
        self.var_set.filter_variants()

        this_consensus = []
        for record in self.var_set:

            this_consensus.append(PHEFilterBase.call_concensus(record))

        self.assertListEqual(expected_consensus, this_consensus)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
