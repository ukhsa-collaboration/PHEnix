"""Prepare reference for different mappers and variant callers."""
#!/usr/bin/env python

import argparse
import logging

from phe.mapping.mapping_factory import available_mappers
from phe.mapping.mapping_factory import factory as map_fac
from phe.variant.variant_factory import available_callers
from phe.variant.variant_factory import factory as var_fac


def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--reference", "-r", required=True, help="Path to reference file to prepare.")

    args.add_argument("--mapper", required=True, help="Available mappers: %s" % available_mappers())
    args.add_argument("--variant", required=True, help="Available variants: %s" % available_callers())

    return args.parse_args()

def main():

    args = get_args()

    mapper = map_fac(args.mapper)
    variant = var_fac(args.variant)

    if not mapper.create_aux_files(args.reference):
        logging.error("Auxiliary files for %s mapper could not be created", args.mapper)

    if not variant.create_aux_files(args.reference):
        logging.error("Auxiliary files for %s variant caller could not be created", args.variant)

    return 0

if __name__ == "__main__":
    exit(main())
