#!/usr/bin/env python
'''
:Date: 7 Apr, 2016
:Author: Alex Jironkin
'''
from distutils.core import setup
import os
import sys

import pip
from pip.req import parse_requirements

def get_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "VERSION")
    version = "N/A"

    if os.path.exists(version_file):
        try:
            with open(version_file) as fp:
                version = fp.next().strip()
                version += "-static"
        except IOError:
            pass
    return version

# At the time of writing there is an open issue on pip > 6.0
#    Where session is required parameter. Breaks backwards compatibility.
if int(pip.__version__.split(".")[0]) >= 6:
    install_reqs = parse_requirements('requirements.txt', session=False)
else:
    install_reqs = parse_requirements('requirements.txt')

install_requires = [str(ir.req) for ir in install_reqs]

setup(name='PHEnix',
      version=get_version(),
      description='Public Health England(UK) SNP calling pipeline tools.',
      author='Public Health England',
      author_email='NGSSBioinformatics@phe.gov.uk',
      url='http://phoenix.readthedocs.org/en/latest/index.html',
      package=['phe',
               'phe.mapping',
               'phe.variant',
               'phe.variant_filters',
               "phe.metadata",
               "phe.annotations"],
      scripts=['scripts/run_snp_pipeline.py',
               "scripts/filter_vcf.py",
               "scripts/prepare_reference.py",
               "scripts/vcf2fasta.py"],
      install_requires=install_requires
     )
