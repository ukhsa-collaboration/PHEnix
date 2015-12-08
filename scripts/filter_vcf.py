#!/usr/bin/env python
'''Simple VCF parser using custom filters.

Created on 6 Oct 2015

@author: alex
'''
import argparse
import logging
import yaml

from phe.variant import VariantSet


def get_args():

    args = argparse.ArgumentParser()

    args.add_argument("--vcf", "-v", required=True, help="VCF file to (re)filter.")

    group = args.add_mutually_exclusive_group()

    group.add_argument("--filters", "-f", help="Filter(s) to apply as key:threshold pairs, separated by comma.")
    group.add_argument("--config", "-c", help="Config with filters in YAML format. E.g.filters:-key:value")

    args.add_argument("--output", "-o", required=True, help="Location for filtered VCF to be written.")

    args.add_argument("--only-good", action="store_true", default=False, help="Write only variants that PASS all filters (default all variants are written).")

    args.add_argument("--debug", action="store_true", default=False, help="Make output more verbose.")

    return args.parse_args()

def load_config(config_path):
    with open(config_path) as fp:
        config = yaml.load(fp)

    return config.get("filters", {})

def main():
    args = get_args()

    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                            level=log_level)

    if args.config is not None:
        args.filters = load_config(args.config)
    elif args.filters is None:
        logging.error("Either --config or --filters needs to be specified.")
        return 1

    var_set = VariantSet(args.vcf, filters=args.filters)

    if args.filters:
        var_set.filter_variants()

    var_set.write_variants(args.output, only_good=args.only_good)



if __name__ == '__main__':
    exit(main())
