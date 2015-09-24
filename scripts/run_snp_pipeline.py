import argparse
import os
import tempfile

import vcf

from phe.mapping.mapping_factory import factory as map_fac
from phe.variant.variant_factory import factory as variant_fac
from phe.variant_filters import filter_vcf, make_filters

def pipeline():
    return 0

def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--workflow", "-w")
    args.add_argument("--input", "-i")

    args.add_argument("-r1")
    args.add_argument("-r2")
    args.add_argument("-r")
    args.add_argument("--outdir", "-o")

    args.add_argument("--mapper", "-m", default="bwa")
    args.add_argument("--variant", "-v", default="gatk")

    return args.parse_args()

def main():
    args = get_args()

    mapper = map_fac(mapper=args.mapper)
    variant = variant_fac(variant=args.variant)

    with tempfile.NamedTemporaryFile(suffix=".sam") as tmp:
        mapper.make_sam(ref=args.r, R1=args.r1, R2=args.r2, out_file=tmp.name)

        os.system("samtools view -bhS %s | samtools sort - %s" % (tmp.name, "test"))
        os.system("samtools index %s" % "test.bam")

#         vcf_file = os.path.abspath("test.vcf")

        vcf_file = variant.make_vcf(ref=args.r, bam="test.bam", out_dir=args.outdir)

        # Example filter options/
        filter_options = {"min_depth": 6, "ad_ratio": 0.89, "gq_score": 1, "mq_score": 31}

        filters = make_filters(config=filter_options)

        good_var, bad_var, template = filter_vcf(vcf_file, filters)

        # Write records to VCF file.
        with open("filtered.vcf", "w") as out_vcf:
            writer = vcf.Writer(out_vcf, template)
            for record in bad_var:
                writer.write_record(record)

    return 0

if __name__ == "__main__":
    exit(main())
