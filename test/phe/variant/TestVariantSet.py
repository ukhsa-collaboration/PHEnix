'''
Created on 3 Dec 2015

@author: alex
'''
import os
import unittest

import vcf

from phe.variant import VariantSet


class TestVariantSet(unittest.TestCase):

    def setUp(self):

        self.config = {"ad_ratio": 0.9, "min_depth": 5, "mq_score": 30}

        base_path = os.path.abspath(os.path.dirname(__file__))
        self.vcf_in = os.path.join(base_path, "..", "variant_filters", "sample.vcf")
        self.vcf_gold_standard = os.path.join(base_path, "sample_gold.vcf")
        self.vcf_snp_only_gold_standard = os.path.join(base_path, "sample_only_good_snps.vcf")
        self.vcf_out = os.path.join(base_path, "sample_out.vcf")

        self.good_vars = [133]
        self.bad_vars = [1, 142, 29144, 31809, 31810, 65032, 65436]

        self.var_set = VariantSet(self.vcf_in, filters=self.config)


    def tearDown(self):
        if os.path.exists(self.vcf_out):
            os.unlink(self.vcf_out)

    def test_filter(self):

        self.var_set.filter_variants()

        good_vars = []
        bad_vars = []

        for var in self.var_set.variants(only_good=False):
            if not var.FILTER:
                good_vars.append(var.POS)
            else:
                bad_vars.append(var.POS)

        good_vars.sort()
        bad_vars.sort()

        self.assertListEqual(self.good_vars, good_vars)
        self.assertListEqual(self.bad_vars, bad_vars)


    def test_write_variants(self):
        self.var_set.filter_variants()

        self.var_set.write_variants(self.vcf_out)

        with open(self.vcf_gold_standard) as fp:
            gold_standard = fp.readlines()

        with open(self.vcf_out) as fp:
            my_file = fp.readlines()


        self.assertListEqual(gold_standard, my_file)

        self.var_set.write_variants(self.vcf_out, only_good=True)

        with open(self.vcf_out) as fp:
            my_file = fp.readlines()
        with open(self.vcf_snp_only_gold_standard) as fp:
            gold_standard = fp.readlines()

        self.assertListEqual(gold_standard, my_file)

    def test_add_metadata(self):
        self.var_set.add_metadata({"key": [{"value": "description"}]})
        self.var_set.write_variants(self.vcf_out)

        vcf_reader = vcf.Reader(filename=self.vcf_out)

        self.assertDictContainsSubset({"key": [{"value": "description"}]}, vcf_reader.metadata)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
