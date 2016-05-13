'''
Created on 3 Dec 2015

@author: alex
'''
from collections import defaultdict
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

        self.good_vars = {"CHR1": [10], "CHR2": [10]}
        self.bad_vars = {"CHR1": [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16], "CHR2": [1, 2, 3, 4, 5, 6, 7, 8, 9]}

        self.var_set = VariantSet(self.vcf_in, filters=self.config)

        for v in self.bad_vars.itervalues():
            v.sort()
        for v in self.good_vars.itervalues():
            v.sort()

    def tearDown(self):
        if os.path.exists(self.vcf_out):
            os.unlink(self.vcf_out)

    def test_filter(self):

        self.var_set.filter_variants()

        good_vars = defaultdict(list)
        bad_vars = defaultdict(list)

        for var in self.var_set.variants(only_good=False):
            if not var.FILTER:
                good_vars[var.CHROM].append(var.POS)
            else:
                bad_vars[var.CHROM].append(var.POS)

        for v in bad_vars.itervalues():
            v.sort()
        for v in good_vars.itervalues():
            v.sort()

        self.assertDictEqual(self.good_vars, good_vars)
        self.assertDictEqual(self.bad_vars, bad_vars)


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
