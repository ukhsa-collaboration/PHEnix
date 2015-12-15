'''
Created on 3 Dec 2015

@author: alex
'''
import os
import unittest
import vcf

from phe.variant_filters.DP4Filter import DP4Filter


class TestADFilter(unittest.TestCase):


    def setUp(self):
        base_path = os.path.abspath(os.path.dirname(__file__))
        self.vcf_in = os.path.join(base_path, "sample.vcf")
        self.filter_threshold = 0.9
        self.parameter = "dp4_ratio"
        self.filter_config = {self.parameter: self.filter_threshold}
        self.filter = DP4Filter(self.filter_config)

        self.bad_positions = [1, 133, 29144, 31809, 31810, 65032, 65436]
        self.good_positions = [142]
        self.na_positions = [133, 29144, 31809, 31810, 65032, 65436]

        self.bad_positions.sort()
        self.na_positions.sort()
        self.good_positions.sort()

    def tearDown(self):
        pass

    def test_call(self):

        reader = vcf.Reader(filename=self.vcf_in)

        bad_positions = []
        good_positions = []
        na_positions = []
        for record in reader:
            result = self.filter(record)

            if result is None:
                good_positions.append(record.POS)
                continue
            elif result is False:
                na_positions.append(record.POS)

            bad_positions.append(record.POS)

        bad_positions.sort()
        good_positions.sort()
        na_positions.sort()

        self.assertListEqual(self.bad_positions, bad_positions)
        self.assertListEqual(self.good_positions, good_positions)
        self.assertListEqual(self.na_positions, na_positions)

    def test_short_desc(self):
        short_desc = "Filter sites by DP4 ratio. (DP4 ratio > %s )" % self.filter_threshold

        self.assertEquals(short_desc, self.filter.short_desc())

    def test_parameter(self):
        self.assertEquals(self.parameter, self.filter.parameter)

    def test_get_config(self):
        self.assertDictEqual(self.filter_config, self.filter.get_config())

    def test_filter_name(self):
        self.assertEquals("%s:%s" % (self.parameter, self.filter_threshold), self.filter.filter_name())

        self.assertEquals("%s:%s" % (self.parameter, self.filter_threshold), str(self.filter))

    def test_bad_config(self):
        with self.assertRaises(Exception):
            DP4Filter({self.parameter: "test"})

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
