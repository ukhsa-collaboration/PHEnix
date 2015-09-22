import argparse
import os
import tempfile

import pysam

from phe.mapping.mapping_factory import factory as map_fac


def pipeline():
    return 0

def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--workflow", "-w")
    args.add_argument("--input", "-i")

    args.add_argument("-r1")
    args.add_argument("-r2")
    args.add_argument("-r")
    args.add_argument("-o")

    args.add_argument("--mapper", "-m", default="bwa")
    args.add_argument("-v", default="gatk")

    return args.parse_args()

def main():
    args = get_args()

    mapper = map_fac(mapper=args.mapper)

    with tempfile.NamedTemporaryFile(suffix=".sam") as tmp:
        mapper.make_sam(ref=args.r, R1=args.r1, R2=args.r2, out_file=tmp.name)

        os.system("samtools view -bhS %s > %s" % (tmp.name, "test.bam"))

    return 0

if __name__ == "__main__":
    exit(main())
