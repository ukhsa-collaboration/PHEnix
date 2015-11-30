'''
Created on 10 Nov 2015

@author: alex
'''
from collections import OrderedDict
import math

import vcf

from phe.annotations import Annotator


class CoverageAnnotator(Annotator):
    '''
    classdocs
    '''

    mean = None
    dev = None
    name = "coverage"

    def __init__(self, config=None):
        '''
        Constructor
        '''
        super(CoverageAnnotator, self).__init__(self.name)

    def annotate(self, vcf_path=None):
        records = []
        dp = []
        reader = vcf.Reader(filename=vcf_path)

        for record in reader:
            records.append(record)

            dp.append(record.INFO.get("DP"))


        self.mean = sum(dp) * 1.0 / len(dp)
        mean_sqr = sum([i ** 2 for i in dp]) * 1.0 / len(dp)
        self.dev = math.sqrt(mean_sqr - self.mean ** 2)

    def get_meta_values(self):
        return OrderedDict({"mean": "%.2f" % self.mean, "dev": "%.2f" % self.dev})
