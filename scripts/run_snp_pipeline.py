#!/usr/bin/env python

from argparse import RawTextHelpFormatter
import argparse
import logging
import os
import sys
import tempfile

import vcf
import yaml

from phe.annotations import available_annotators, make_annotators
from phe.mapping.mapping_factory import factory as map_fac, available_mappers
from phe.variant import VariantSet
from phe.variant.variant_factory import factory as variant_fac, \
    available_callers
from phe.variant_filters import available_filters, str_to_filters, make_filters


def pipeline():
    return 0

desc = '''Run the snp pipeline with specified mapper, variant caller and some filters.
Available mappers: %s

Available variant callers: %s

Available filters: %s

Available annotators: %s''' % (available_mappers(), available_callers(), available_filters(), available_annotators())

def get_args():
    args = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)

    args.add_argument("--workflow", "-w")
    args.add_argument("--input", "-i")

    args.add_argument("-r1", help="R1/Forward read in Fastq format.")
    args.add_argument("-r2", help="R2/Reverse read in Fastq format.")
    args.add_argument("-r", help="Rerefence to use for mapping.")
    args.add_argument("--sample-name", default="test_sample", help="Name of the sample for mapper to include as read groups.")
    args.add_argument("--outdir", "-o")

    args.add_argument("--config", "-c")

    args.add_argument("--mapper", "-m", default="bwa", help="Available mappers: %s" % available_mappers())
    args.add_argument("--mapper-options", help="Custom maper options (advanced)")
    args.add_argument("--bam")
    args.add_argument("--variant", "-v", default="gatk", help="Available variant callers: %s" % available_callers())
    args.add_argument("--variant-options", help="Custom variant options (advanced)")
    args.add_argument("--vcf")
    args.add_argument("--filters", type=str, help="Filters to be applied to the VCF in key:value pairs, separated by comma (,). Available_filters: %s" % available_filters())

    args.add_argument("--annotators", nargs="+", help="List of annotators to run before filters. Available: %s" % available_annotators())

    return args.parse_args()

def load_config(args):

    with open(args.config) as fp:
        config = yaml.load(fp)

    args.mapper = config.get("mapper")
    args.mapper_options = config.get("mapper-options")

    args.variant = config.get("variant")
    args.variant_options = config.get("variant-options")

    args.filters = config.get("filters")


def main():

    logging.basicConfig(level=logging.DEBUG,)

    logging.info("Initialising data matrix.")
    args = get_args()

    if args.outdir is None:
        sys.stdout.write("Please provide output directory.")
        return -1
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # If config is specified, then load data from that.
    if args.config:
        load_config(args)

    mapper = map_fac(mapper=args.mapper, custom_options=args.mapper_options)

    variant = None
    if args.variant:
        variant = variant_fac(variant=args.variant, custom_options=args.variant_options)

    if args.annotators:
        args.annotators = make_annotators(args.annotators)

    if args.filters:
        try:
            if isinstance(args.filters, str):
                args.filters = str_to_filters(args.filters)
            elif isinstance(args.filters, dict):
                args.filters = make_filters(args.filters)
            else:
                logging.warn("Unknown filters specified: %s", args.filters)
        except Exception:
            logging.error("Failed to recognise and create filters.")
            return 3

    logging.info("Mapping data file.")
    if args.bam is not None:
        bam_file = args.bam
    elif args.vcf is None:
        bam_file = os.path.join(args.outdir, "%s.bam" % args.sample_name)
        success = mapper.make_bam(ref=args.r, R1=args.r1, R2=args.r2, out_file=bam_file, sample_name=args.sample_name)

        if not success:
            logging.warn("Could not map reads to the reference. Aborting.")
            return 1

    logging.info("Creating digitised variants.")
    if args.vcf:
        vcf_file = args.vcf
    else:
        vcf_file = os.path.join(args.outdir, "%s.vcf" % args.sample_name)

        if variant and not variant.make_vcf(ref=args.r, bam=bam_file, vcf_file=vcf_file):
            logging.error("VCF was not created.")
            return 2

    annotators_metadata = []
    if args.annotators:
        logging.info("Annotating")
        for annotator in args.annotators:
            # TODO: This iterates over the VCF for each annotator. Not good.
            annotator.annotate(vcf_path=vcf_file)

            meta = annotator.get_meta()

            if meta:
                annotators_metadata.append(meta)

    if args.filters:
        logging.info("Applying filters: %s", args.filters)
        var_set = VariantSet(vcf_file, filters=args.filters)

        var_set.add_metadata(mapper.get_meta())
        var_set.add_metadata(variant.get_meta())

        for annotator_md in annotators_metadata:
            var_set.add_metadata(annotator_md)

        var_set.filter_variants()

#         var_set.write_variants("filtered.vcf", only_snps=True, only_good=True)
#
#         var_set.write_variants("filtered.all.vcf")
#
#         var_set._write_bad_variants("filtered.bad.vcf")

        final_vcf = os.path.join(args.outdir, "%s.filtered.vcf" % args.sample_name)
        var_set.serialise(final_vcf)

    return 0

if __name__ == "__main__":
    exit(main())
