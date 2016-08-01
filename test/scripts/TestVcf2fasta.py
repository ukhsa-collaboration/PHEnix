'''
Created on 1 Aug 2016

@author: ulf
'''
import sys
import os
import unittest
import tempfile
import shutil

import vcf

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts')))
import vcf2fasta

class TestVcf2fasta(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.abspath(os.path.dirname(__file__))
        self.fasta_gold_std = os.path.join(self.base_dir, "data", "vcf2fasta_gold.fasta")
        self.outfolder = os.path.join(tempfile.mkdtemp())
        self.args = {"input": [os.path.join(self.base_dir, "data", x) for x in ["test1.vcf", "test2.vcf"]],
                     "out": os.path.join(self.outfolder, "result.fasta"),
                     "tmp": None,
                     "exclude": None,
                     "include": None,
                     "directory": None,
                     "with_stats": None,
                     "column_Ns": None,
                     "sample_Ns": None,
                     "column_gaps": None,
                     "sample_gaps": None,
                     "reference": None}

    # ----------------------------------------------------------------------------------------------

    def tearDown(self):
        shutil.rmtree(self.outfolder)

    # ----------------------------------------------------------------------------------------------

    def test_vcf2fasta(self):

        vcf2fasta.main(self.args)

        with open(self.fasta_gold_std, 'rb') as fp:
            gold_standard = fp.readlines()

        with open(self.args['out'], 'rb') as fp:
            my_file = fp.readlines()

        self.assertListEqual(gold_standard, my_file)

# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
