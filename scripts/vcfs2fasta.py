'''
Created on 5 Oct 2015

@author: alex
'''
import argparse
import os

import vcf

from phe.variant_filters import PHEFilterBase, make_filters


def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--input", "-i", type=str, nargs='+')
    args.add_argument("--out", "-o")


    return args.parse_args()

def main():

    args = get_args()

    all_pos = dict()
    ref = None

    for vcf_in in args.input:
        reader = vcf.Reader(filename=vcf_in)

        for record in reader:
            if ref is None:
                ref = record.CHROM
            else:
                if ref != record.CHROM:
                    print "Chromosomes don't match."
                    return 1
            if record.FILTER == "PASS" or not record.FILTER:
                if record.is_snp:
                    if record.POS in all_pos and all_pos[record.POS] != record.REF:
                        print "SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s" % record.POS
                        return 2

                    all_pos[record.POS] = record.REF


    all_data = {}

    for vcf_in in args.input:
        reader = vcf.Reader(filename=vcf_in)

        sample_seq = ""
        sample_name, _ = os.path.splitext(vcf_in)

        all_data[sample_name] = { pos: all_pos[pos] for pos in all_pos.keys() }

        for record in reader:
            if record.POS in all_pos.keys():
                if record.FILTER == "PASS" or not record.FILTER:
                    if record.is_snp:
                        if len(record.ALT) > 1:
                            print "POS %s passed filters but has multiple alleles. Inserting N"
                            all_data[sample_name][record.POS] = "N"
                        else:
                            all_data[sample_name][record.POS] = record.ALT[0].sequence
                else:
                    filters = {}
                    for filter_id in record.FILTER:
                        filters.update(PHEFilterBase.decode(filter_id))

                    if filters:
                        filters = make_filters(config=filters)

                    is_gap = [True for f in filters if f.is_gap()]
                    if len(is_gap) > 0:
                        all_data[sample_name][record.POS] = "-"
                    else:
                        all_data[sample_name][record.POS] = "N"


    with open(args.out, "w") as fp:
        for sample in all_data:
            sample_seq = ""
            for pos in all_pos.keys():
                sample_seq += all_data[sample][pos]

            if len(sample_seq) == len(all_pos.keys()):
                fp.write(">%s\n%s\n" % (sample_name, sample_seq))
    return 0

if __name__ == '__main__':
    exit(main())
