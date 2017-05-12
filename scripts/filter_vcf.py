'''Simple VCF parser using custom filters.

:Date: 6 October, 2015
:Author: alex
'''
import argparse
from collections import OrderedDict
import logging
import os

import yaml

from phe.variant import VariantSet


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

def get_desc():
    return "Filter the VCF using provided filters."

def get_args():
    args = argparse.ArgumentParser(description=get_desc())

    args.add_argument("--vcf", "-v", required=True, help="VCF file to (re)filter.")

    group = args.add_mutually_exclusive_group()

    group.add_argument("--filters",
                       "-f",
                       help="Filter(s) to apply as key:threshold pairs, separated by comma. \
                            Recommendations: GATK: mq_score:30,min_depth:10,ad_ratio:0.9 \
                                             Mpileup: mq_score:30,min_depth:10,dp4_ratio:0.9")
    group.add_argument("--config", "-c", help="Config with filters in YAML format. E.g.filters:-key:value")

    args.add_argument("--output", "-o", required=True, help="Location for filtered VCF to be written.")

    args.add_argument("--reference", "-r", help="mpileup version <= 1.3 do not output all positions. This is required to fix rfrence base in VCF.")

    args.add_argument("--only-good", action="store_true", default=False, help="Write only variants that PASS all filters (default all variants are written).")

    return args

def load_config(config_path):
    with open(config_path) as fp:
        config = yaml.load(fp)

    return config.get("filters", {})

def main(args):

    if args.get("version") is None:
        args["version"] = get_version()

    if args["config"] is not None:
        args["filters"] = load_config(args["config"])
    elif args["filters"] is None and not args["only_good"]:
        logging.error("Either --config or --filters needs to be specified.")
        return 1

    var_set = VariantSet(args["vcf"], filters=args["filters"], reference=args["reference"])

    if args.get("version") is not None:
        var_set.add_metadata(OrderedDict({"PHEnix-Version": (args["version"],)}))

    var_set.filter_variants(out_vcf=args["output"], only_good=args["only_good"])

    logging.info("Finished filtering")
    return 0

if __name__ == '__main__':
    exit(main(vars(get_args().parse_args())))
