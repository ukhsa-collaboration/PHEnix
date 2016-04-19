'''Simple VCF parser using custom filters.

:Date: 6 October, 2015
:Author: alex
'''
import argparse
from collections import OrderedDict
import logging

import yaml

from phe.variant import VariantSet
import phenix_versioneer


def get_desc():
    return "Filter the VCF using provided filters."

def get_args():
    args = argparse.ArgumentParser(description=get_desc())

    args.add_argument("--vcf", "-v", required=True, help="VCF file to (re)filter.")

    group = args.add_mutually_exclusive_group()

    group.add_argument("--filters", "-f", help="Filter(s) to apply as key:threshold pairs, separated by comma.")
    group.add_argument("--config", "-c", help="Config with filters in YAML format. E.g.filters:-key:value")

    args.add_argument("--output", "-o", required=True, help="Location for filtered VCF to be written.")

    args.add_argument("--only-good", action="store_true", default=False, help="Write only variants that PASS all filters (default all variants are written).")

    return args

def load_config(config_path):
    with open(config_path) as fp:
        config = yaml.load(fp)

    return config.get("filters", {})

def main(args):

    if args.get("version") is None:
        args["version"] = phenix_versioneer.get_version()

    if args["config"] is not None:
        args["filters"] = load_config(args["config"])
    elif args["filters"] is None and not args["only_good"]:
        logging.error("Either --config or --filters needs to be specified.")
        return 1

    var_set = VariantSet(args["vcf"], filters=args["filters"])

    if args.get("version") is not None:
        var_set.add_metadata(OrderedDict({"PHEnix-Version": (args["version"],)}))

    if args["filters"]:
        var_set.filter_variants()

    var_set.write_variants(args["output"], only_good=args["only_good"])

    logging.info("Finished filtering")
    return 0

if __name__ == '__main__':
    exit(main(vars(get_args().parse_args())))
