'''
:Date: 10 Nov, 2015
:Author: Alex Jironkin
'''
from collections import OrderedDict
import math

import vcf

from phe.annotations import Annotator


class CoverageAnnotator(Annotator):
    '''
    classdocs
    '''

    name = "coverage"

    def __init__(self, config=None):
        '''Constructor'''
        super(CoverageAnnotator, self).__init__(self.name)

        self.mean = 0
        self._mean_sqr = 0
        self.dev = 0

    def annotate(self, vcf_path=None):
        reader = vcf.Reader(filename=vcf_path)
        total = 0
        for record in reader:

            self.mean += record.INFO.get("DP", 0)
            self._mean_sqr += record.INFO.get("DP", 0) ** 2
            total += 1


        self.mean = self.mean * 1.0 / total
        self._mean_sqr = self._mean_sqr * 1.0 / total
        self.dev = math.sqrt(self._mean_sqr - self.mean ** 2)

    def get_meta_values(self):
        return OrderedDict({"mean": "%.2f" % self.mean, "dev": "%.2f" % self.dev})
