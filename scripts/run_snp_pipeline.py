from argparse import RawTextHelpFormatter
import argparse
import logging
import os
import tempfile

import vcf

from phe.mapping.mapping_factory import factory as map_fac, available_mappers
from phe.variant import VariantSet
from phe.variant.variant_factory import factory as variant_fac, \
    available_callers
from phe.variant_filters import available_filters


def pipeline():
    return 0

desc = '''Run the snp pipeline with specified mapper, variant caller and some filters.
Available mappers: %s

Available variant callers: %s

Available filters: %s''' % (available_mappers(), available_callers(), available_filters())

def get_args():
    args = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)

    args.add_argument("--workflow", "-w")
    args.add_argument("--input", "-i")

    args.add_argument("-r1")
    args.add_argument("-r2")
    args.add_argument("-r")
    args.add_argument("--sample_name", default="test_sample", help="Name of the sample for mapper to include as read groups.")
    args.add_argument("--outdir", "-o")

    args.add_argument("--mapper", "-m", default="bwa", help="Available mappers: %s" % available_mappers())
    args.add_argument("--mapper-options", help="Custom maper options (advanced)")
    args.add_argument("--variant", "-v", default="gatk", help="Available variant callers: %s" % available_callers())
    args.add_argument("--variant-options", help="Custom variant options (advanced)")
    args.add_argument("--filters", type=str, help="Filters to be applied to the VCF in key:value pairs, separated by comma (,). Available_filters: %s" % available_filters())

    return args.parse_args()

def main():

    logging.basicConfig(level=logging.DEBUG,)

    logging.info("Initialising data matrix.")
    args = get_args()

    mapper = map_fac(mapper=args.mapper, custom_options=args.mapper_options)
    variant = variant_fac(variant=args.variant, custom_options=args.variant_options)

    logging.info("Computing registry information.")
    if args.filters is not None:
        filters = {}
        for kv_pair in args.filters.split(","):
            pair = kv_pair.split(":")
            assert len(pair) == 2, "Filters should be separated by ':' %s" % kv_pair

            # We don't care about casting them to correct type because Filters
            #    will do it for us.
            filters[pair[0]] = pair[1]

        args.filters = filters


    with tempfile.NamedTemporaryFile(suffix=".bam") as tmp:
        logging.info("Mapping data file.")
        success = mapper.make_bam(ref=args.r, R1=args.r1, R2=args.r2, out_file="%s.bam" % args.sample_name, sample_name=args.sample_name)

        if not success:
            logging.warn("Could not map reads to the reference. Aborting.")
            return 1

#         vcf_file = os.path.abspath("all_variants.vcf")
        logging.info("Creating digitised variants.")

        vcf_file = variant.make_vcf(ref=args.r, bam="%s.bam" % args.sample_name, out_dir=args.outdir)

        if args.filters:
            logging.info("Applying filters: %s", args.filters)
            var_set = VariantSet(vcf_file, args.filters)

            var_set.make_variant_set()

            var_set.write_variants("filtered.vcf", only_snps=True)

            var_set._write_bad_variants("filtered.bad.vcf")

#             var_set.serialise("var_set.pkl")

    return 0

if __name__ == "__main__":
    exit(main())
