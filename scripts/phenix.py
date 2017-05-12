#!/usr/bin/env python
'''
Single script for running other scripts in the project.

:Date: 14 April, 2016
:Author: Public Health England
'''

from argparse import RawTextHelpFormatter
import argparse
import logging
import os

import filter_vcf
import prepare_reference
import run_snp_pipeline
import vcf2fasta
import vcf2distancematrix
import vcf2json


def get_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "VERSION")
    version = "N/A"

    if os.path.exists(version_file):
        try:
            with open(version_file) as fp:
                version = fp.next().strip()
        except IOError:
            pass
    return version

def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--debug", action="store_true", help="More verbose logging (default: turned off).")

    args.add_argument("--version", action="version", version=get_version())

    subparsers = args.add_subparsers(dest='cmd')

    # NOTES: add_help=False is also required, otherwise args conflict
    #    error will be shown.

    subparsers.add_parser("run_snp_pipeline",
                          description=run_snp_pipeline.get_desc(),
                          help="Run SNP pipeline.",
                          parents=[run_snp_pipeline.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("filter_vcf",
                          description=filter_vcf.get_desc(),
                          help="Filter a VCF.",
                          parents=[filter_vcf.get_args()],
                          add_help=False)

    subparsers.add_parser("prepare_reference",
                          description=prepare_reference.get_desc(),
                          help="Create aux files for reference.",
                          parents=[prepare_reference.get_args()],
                          add_help=False)

    subparsers.add_parser("vcf2fasta",
                          description=vcf2fasta.get_desc(),
                          help="Convert VCFs to FASTA.",
                          parents=[vcf2fasta.get_args()],
                          add_help=False)

    subparsers.add_parser("vcf2distancematrix",
                          description=vcf2distancematrix.get_desc(),
                          help="Convert VCFs to a distance matrix.",
                          parents=[vcf2distancematrix.get_args()],
                          add_help=False)

    subparsers.add_parser("vcf2json",
                          description=vcf2json.get_desc(),
                          help="Convert VCFs to a JSON file containing variants and ignored positions as arrays of positions relative to reference chromosomes.",
                          parents=[vcf2json.get_args()],
                          add_help=False)

    return args

def main():
    version = get_version()

    args = vars(get_args().parse_args())

    log_level = logging.DEBUG if args["debug"] else logging.INFO
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                            level=log_level)

    logging.info("Version: %s", version)

    args["version"] = version

    if args["cmd"] == "run_snp_pipeline":
        return run_snp_pipeline.main(args)
    elif args["cmd"] == "filter_vcf":
        return filter_vcf.main(args)
    elif args["cmd"] == "prepare_reference":
        return prepare_reference.main(args)
    elif args["cmd"] == "vcf2fasta":
        return vcf2fasta.main(args)
    elif args["cmd"] == "vcf2distancematrix":
        return vcf2distancematrix.main(args)
    elif args["cmd"] == "vcf2json":
        return vcf2json.main(args)
    return 1

if __name__ == '__main__':
    exit(main())
