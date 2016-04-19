"""Prepare reference for different mappers and variant callers.

:Date: 6 October, 2015
:Author: Alex Jironkin
"""

import argparse
import logging

from phe.mapping.mapping_factory import available_mappers
from phe.mapping.mapping_factory import factory as map_fac
from phe.variant.variant_factory import available_callers
from phe.variant.variant_factory import factory as var_fac

def get_desc():
    return "Prepare reference for SNP pipeline by generating required aux files."

def get_args():
    args = argparse.ArgumentParser(description=get_desc())

    args.add_argument("--reference", "-r", required=True, help="Path to reference file to prepare.")

    args.add_argument("--mapper", help="Available mappers: %s" % available_mappers())
    args.add_argument("--variant", help="Available variants: %s" % available_callers())

    return args

def main(args):

    result = 0

    logging.info("Creating auxilliary files for %s", args["reference"])

    if args["mapper"]:
        mapper = map_fac(args["mapper"])

        if mapper is None or not mapper.create_aux_files(args["reference"]):
            logging.error("Auxiliary files for %s mapper could not be created", args["mapper"])
            result += 1

    if args["variant"]:
        variant = var_fac(args["variant"])
        if variant is None or not variant.create_aux_files(args["reference"]):
            logging.error("Auxiliary files for %s variant caller could not be created", args["variant"])
            result += 1

    if args["mapper"] is None and args["variant"] is None:
        logging.info("Nothing specified. Please provider either --mapper or --variant.")

    logging.info("Finished creating auxilliary files.")

    return result

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
