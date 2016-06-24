'''
Created on 22 Apr 2016

@author: alex
'''
import logging
from math import ceil
import os
import shutil
import sys
import unittest

import vcf


sys.path.append("../../scripts")

class RunSnpPipeline(unittest.TestCase):

    def check_common_files(self):

        base_dir = self.args["outdir"]
        self.assertTrue(os.path.exists(base_dir))

        self.assertTrue(os.path.exists(os.path.join(base_dir, "test_sample.vcf")), "Missing RAW vcf: %s" % base_dir)
        self.assertTrue(os.path.exists(os.path.join(base_dir, "test_sample.filtered.vcf")), "Missing FILTERED vcf: %s" % base_dir)
        self.assertTrue(os.path.exists(os.path.join(base_dir, "test_sample.bam")), "Missing BAM file: %s" % base_dir)
        self.assertTrue(os.path.exists(os.path.join(base_dir, "test_sample.bam.bai")), "Missing BAM INDEX file: %s" % base_dir)

    def check_variants(self, vcf_file):
        reader = vcf.Reader(filename=vcf_file)

        missing_snps = {}
        missing_bad = {}
        for r in reader:
            if not r.FILTER:
                try:
                    if r.is_snp:
                        self.good_vars.remove(r.POS)
                except KeyError:
                    for f in r.FILTER:
                        if f not in missing_snps:
                            missing_snps[f] = 0
                        missing_snps[f] += 1
            else:
                try:
                    self.failed_vars.remove(r.POS)
                except KeyError:
                    for f in r.FILTER:
                        if f not in missing_bad:
                            missing_bad[f] = 0
                        missing_bad[f] += 1



        self.assertEqual(missing_snps, {})

        self.assertLessEqual(len(self.good_vars), self.good_threshold,
                             ">=%.2f different good positions: [%s]" % (self.threshold, ",".join(str(i) for i in self.failed_vars)))

        self.assertLessEqual(len(self.failed_vars),
                             self.failed_threshold,
                             ">=%.2f failed positions: [%s]" % (self.threshold, ",".join(str(i) for i in self.failed_vars)))

    def setUp(self):

        logging.basicConfig(level=logging.DEBUG)

        self.base_dir = os.path.abspath(os.path.dirname(__file__))
        self.args = {"r1": os.path.join(self.base_dir, "data", "R1.fastq.gz"),
                "r2": os.path.join(self.base_dir, "data", "R2.fastq.gz"),
                "reference": os.path.join(self.base_dir, "data", "reference.fa"),
                "workflow": None,
                "input": None,
                "outdir": "results",
                "bam": None,
                "vcf": None,
                "sample_name": "test_sample",
                "keep_temp": True
                }

        if not os.path.exists(self.args["outdir"]):
            os.makedirs(self.args["outdir"])

        self.good_vars = set()
        with open(os.path.join(self.base_dir, "data", "good_vars.csv")) as fp:
            for line in fp:
                self.good_vars.add(int(line.strip().split(",")[0]))

        self.failed_vars = set()
        with open(os.path.join(self.base_dir, "data", "failed_vars.csv")) as fp:
            for line in fp:
                self.failed_vars.add(int(line.strip().split(",")[0]))

        self.threshold = 0.05

        self.good_threshold = int(ceil(len(self.good_vars) * self.threshold))
        self.failed_threshold = int(ceil(len(self.failed_vars) * self.threshold))

    def tearDown(self):
        pass
#         if os.path.exists(self.args["outdir"]):
#             shutil.rmtree(self.args["outdir"])

    def testEndToEndBwaGatk(self):

        self.args["outdir"] = os.path.join(self.args["outdir"], "bwa_gatk")
        self.args["config"] = os.path.join(self.base_dir, "data", "bwa_gatk.yml")

        self.check_common_files()

        test_vcf = os.path.join(self.args["outdir"], "test_sample.filtered.vcf")

#         self.check_variants(test_vcf)

    def testEndToEndBowtieGatk(self):

        self.args["outdir"] = os.path.join(self.args["outdir"], "bowtie_gatk")
        self.args["config"] = os.path.join(self.base_dir, "data", "bowtie_gatk.yml")

        self.check_common_files()

        test_vcf = os.path.join(self.args["outdir"], "test_sample.filtered.vcf")

#         self.check_variants(test_vcf)

    def testEndToEndBwaMPileup(self):

        self.args["outdir"] = os.path.join(self.args["outdir"], "bwa_mpileup")

        self.args["config"] = os.path.join(self.base_dir, "data", "bwa_mpileup.yml")

        self.check_common_files()

        test_vcf = os.path.join(self.args["outdir"], "test_sample.filtered.vcf")

#         self.check_variants(test_vcf)

    def testEndToEndBowtieMPileup(self):

        self.args["outdir"] = os.path.join(self.args["outdir"], "bowtie_mpileup")

        self.args["config"] = os.path.join(self.base_dir, "data", "bowtie_mpileup.yml")

        self.check_common_files()

        test_vcf = os.path.join(self.args["outdir"], "test_sample.filtered.vcf")

#         self.check_variants(test_vcf)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
