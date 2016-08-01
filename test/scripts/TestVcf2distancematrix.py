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
import vcf2distancematrix

class TestVcf2distancematrix(unittest.TestCase):

    def setUp(self):
        self.base_dir = os.path.abspath(os.path.dirname(__file__))
        self.tsv_gold_std = os.path.join(self.base_dir, "data", "vcf2distancematrix_gold.tsv")
        self.tre_gold_std = os.path.join(self.base_dir, "data", "vcf2distancematrix_gold.tree")
        self.outfolder = os.path.join(tempfile.mkdtemp())
        self.args = {"input": [os.path.join(self.base_dir, "data", x) for x in ["test1.vcf", "test2.vcf"]],
                     "out": os.path.join(self.outfolder, "result.tsv"),
                     "tree": os.path.join(self.outfolder, "result.tre"),
                     "deletion": None,
                     "format": 'tsv',
                     "substitution": 'number_of_differences',
                     "exclude": None,
                     "include": None,
                     "directory": None,
                     "winsize": 1000,
                     "refgenome": os.path.join(self.base_dir, "data", 'Neisseria_gonorrhoeae_NCCP11945.fasta'),
                     "threads": 1,
                     "with_stats": False,
                     "remove_recombination": True}

    # ----------------------------------------------------------------------------------------------

    def tearDown(self):
        shutil.rmtree(self.outfolder)

    # ----------------------------------------------------------------------------------------------

    def test_vcf2fasta(self):

        vcf2distancematrix.main(self.args)

        with open(self.tsv_gold_std, 'rb') as fp:
            gold_standard = fp.readlines()

        with open(self.args['out'], 'rb') as fp:
            my_file = fp.readlines()

        self.assertListEqual(gold_standard, my_file)


        with open(self.tre_gold_std, 'rb') as fp:
            gold_standard = fp.readlines()

        with open(self.args['tree'], 'rb') as fp:
            my_file = fp.readlines()

        self.assertListEqual(gold_standard, my_file)


# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
