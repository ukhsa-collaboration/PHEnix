'''
Merge SNP data from multiple VCF files into a single fasta file.

Created on 5 Oct 2015

@author: alex
'''
import argparse
import glob
import os

from bintrees import FastRBTree
import vcf

from phe.variant_filters import PHEFilterBase, make_filters


def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--directory", "-d", help="Path to the directory with .vcf files.")
    args.add_argument("--input", "-i", type=str, nargs='+', help="List of VCF files to process.")
    args.add_argument("--out", "-o", required=True, help="Path to the output FASTA file.")


    return args.parse_args()

def main():
    """
    Process VCF files and merge them into a single fasta file.
    """

    args = get_args()


    all_pos = dict()
    contigs = list()
    ref_seq = dict()

    # All positions available for analysis.
    avail_pos = dict()
    # Stats about each position in each chromosome.
    pos_stats = dict()
    # Cached version of the data.
    vcf_data = dict()

    if args.directory is not None and args.input is None:
        args.input = glob.glob(os.path.join(args.directory, "*.vcf"))

    # First pass to get the references and the positions to be analysed.
    for vcf_in in args.input:
        vcf_data[vcf_in] = list()
        reader = vcf.Reader(filename=vcf_in)

        for record in reader:
            vcf_data[vcf_in].append(record)

            if record.CHROM not in contigs:
                contigs.append(record.CHROM)
                avail_pos[record.CHROM] = FastRBTree()

            if record.FILTER == "PASS" or not record.FILTER:
                if record.is_snp:
                    if record.POS in all_pos and all_pos[record.POS] != record.REF:
                        print "SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s" % record.POS
                        return 2

                    if record.CHROM not in pos_stats:
                        pos_stats[record.CHROM] = {}

                    if record.CHROM not in all_pos:
                        all_pos[record.CHROM] = {}

                    if record.CHROM not in ref_seq:
                        ref_seq[record.CHROM] = {}

                    avail_pos[record.CHROM].insert(record.POS, str(record.REF))

                    pos_stats[record.CHROM][record.POS] = {"N":0, "-": 0}
                    all_pos[record.CHROM][record.POS] = record.REF
                    ref_seq[record.CHROM][record.POS] = record.REF


    all_data = { contig: {} for contig in contigs}
    samples = []

    # Instead of loading filters for each position, cache them.
    cached_filters = {}

    for vcf_in in args.input:

        sample_seq = ""
        sample_name, _ = os.path.splitext(os.path.basename(vcf_in))
        samples.append(sample_name)

        # Initialise the data for this sample to be REF positions.
        for contig in contigs:
            all_data[contig][sample_name] = { pos: avail_pos[contig][pos] for pos in avail_pos[contig] }

        for record in vcf_data[vcf_in]:
            # Array of filters that have been applied.
            filters = []

            # If position is our available position.
            if avail_pos[record.CHROM].get(record.POS, False):
                if record.FILTER == "PASS" or not record.FILTER:
                    if record.is_snp:
                        if len(record.ALT) > 1:
                            print "POS %s passed filters but has multiple alleles. Inserting N"
                            all_data[record.CHROM][sample_name][record.POS] = "N"
                        else:
                            all_data[record.CHROM][sample_name][record.POS] = record.ALT[0].sequence
                else:
                    for filter_id in record.FILTER:
                        this_filter = PHEFilterBase.decode(filter_id)

                        if tuple(this_filter.items()) not in cached_filters:
                            cached_filters[tuple(this_filter.items())] = make_filters(config=this_filter)

                        filters += cached_filters[tuple(this_filter.items())]

                    # Currently we are only using first filter to call consensus.
                    extended_code = filters[0].call_concensus(record)

                    # Calculate the stats
                    if extended_code == "N":
                        pos_stats[record.CHROM][record.POS]["N"] += 1
                    elif extended_code == "-":
                        pos_stats[record.CHROM][record.POS]["-"] += 1

                    # Save the extended code of the SNP.
                    all_data[record.CHROM][sample_name][record.POS] = extended_code

    # Output the data to the fasta file.
    # The data is already aligned so simply output it.
    discarded = 0
    with open(args.out, "w") as fp:
        for sample in samples:
            sample_seq = ""
            for contig in contigs:
                for pos in all_pos[contig]:
                    if float(pos_stats[contig][pos]["N"]) / len(all_data[contig]) < 0.1 and \
                        float(pos_stats[contig][pos]["-"]) / len(all_data[contig]) < 0.1:
                        sample_seq += all_data[contig][sample][pos]
                    else:
                        discarded += 1

            fp.write(">%s\n%s\n" % (sample, sample_seq))
        # Do the same for reference data.
        ref_snps = ""
        for contig in contigs:
            for pos in ref_seq[contig]:
                if float(pos_stats[contig][pos]["N"]) / len(all_data[contig]) < 0.1 and \
                        float(pos_stats[contig][pos]["-"]) / len(all_data[contig]) < 0.1:
                    ref_snps += str(ref_seq[contig][pos])
        fp.write(">reference\n%s\n" % ref_snps)

    print("Discarded total of %i from %i for poor quality" % (float(discarded) / len(all_data), len(pos_stats)))
    return 0

if __name__ == '__main__':
    import time

#     with PyCallGraph(output=graphviz):
#     T0 = time.time()
    r = main()
#     T1 = time.time()

#     print "Time taken: %i" % (T1 - T0)
    exit(r)
