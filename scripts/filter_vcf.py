#!/usr/bin/env python
'''Simple VCF parser using custom filters.

Created on 6 Oct 2015

@author: alex
'''
import argparse

from phe.variant import VariantSet


def get_args():

    args = argparse.ArgumentParser()

    args.add_argument("--vcf", "-v", required=True, help="VCF file to (re)filter.")
    args.add_argument("--filters", "-f", required=True, help="Filter(s) to apply as key:threshold pairs, separated by comma.")
    args.add_argument("--output", "-o", required=True, help="Location for filtered VCF to be written.")

    return args.parse_args()


def main():
    args = get_args()

    var_set = VariantSet(args.vcf, filters=args.filters)

    var_set.filter_variants()

    var_set.serialise(args.output)

if __name__ == '__main__':
    exit(main())
