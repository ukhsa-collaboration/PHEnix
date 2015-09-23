import argparse
import os
import tempfile

# import pysam

from phe.mapping.mapping_factory import factory as map_fac
from phe.variant.variant_factory import factory as variant_fac


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

        variant.make_vcf(ref=args.r, bam="test.bam", out_dir=args.outdir)

    return 0

if __name__ == "__main__":
    exit(main())
