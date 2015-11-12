#!/usr/bin/env python
'''
Merge SNP data from multiple VCF files into a single fasta file.

Created on 5 Oct 2015

@author: alex
'''
import argparse
from collections import OrderedDict
import glob
import itertools
import logging
import os

from Bio import SeqIO
from bintrees import FastRBTree
from matplotlib import pyplot as plt
import numpy
import vcf

from phe.variant_filters import IUPAC_CODES


def plot_stats(pos_stats, plots_dir="plots"):
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    for contig in pos_stats:
        x = numpy.array(pos_stats[contig].keys())
        y = numpy.array([ pos_stats[contig][pos]["mut"] for pos in pos_stats[contig]])

        plt.figure()
        plt.plot(x, y, 'ro')

        plt.savefig(os.path.join(plots_dir, "%s.png" % contig))

def get_mixture(record, threshold):
    mixtures = {}
    try:
        if len(record.samples[0].data.AD) > 1:

            total_depth = sum(record.samples[0].data.AD)
            # Go over all combinations of touples.
            for comb in itertools.combinations(range(0, len(record.samples[0].data.AD)), 2):
                i = comb[0]
                j = comb[1]

                alleles = list()

                if 0 in comb:
                    alleles.append(str(record.REF))

                if i != 0:
                    alleles.append(str(record.ALT[i - 1]))
                    mixture = record.samples[0].data.AD[i]
                if j != 0:
                    alleles.append(str(record.ALT[j - 1]))
                    mixture = record.samples[0].data.AD[j]

                ratio = float(mixture) / total_depth
                if ratio == 1.0:
                    logging.debug("This is only designed for mixtures! %s %s %s %s", record, ratio, record.samples[0].data.AD, record.FILTER)

                    if ratio not in mixtures:
                        mixtures[ratio] = []
                    mixtures[ratio].append(alleles.pop())

                elif ratio >= threshold:
                    try:
                        code = IUPAC_CODES[frozenset(alleles)]
                        if ratio not in mixtures:
                            mixtures[ratio] = []
                            mixtures[ratio].append(code)
                    except KeyError:
                        logging.warn("Could not retrieve IUPAC code for %s from %s", alleles, record)
    except AttributeError:
        mixtures = {}

    return mixtures

def print_stats(stats, pos_stats, total_vars):
    for contig in stats:
        for sample, info in stats[contig].items():
            print "%s,%i,%i" % (sample, len(info.get("n_pos", [])), total_vars)

    for contig in stats:
        for pos, info in pos_stats[contig].iteritems():
            print "%s,%i,%i,%i,%i" % (contig, pos, info.get("N", "NA"), info.get("-", "NA"), info.get("mut", "NA"))


def get_args():
    args = argparse.ArgumentParser(description="Combine multiple VCFs into a single FASTA file.")

    group = args.add_mutually_exclusive_group(required=True)
    group.add_argument("--directory", "-d", help="Path to the directory with .vcf files.")
    group.add_argument("--input", "-i", type=str, nargs='+', help="List of VCF files to process.")

    args.add_argument("--out", "-o", required=True, help="Path to the output FASTA file.")

    args.add_argument("--with-mixtures", type=float, help="Specify this option with a threshold to output mixtures above this threshold.")

    args.add_argument("--column-Ns", type=float, help="Keeps columns with fraction of Ns above specified threshold.")

    args.add_argument("--sample-Ns", type=float, help="Keeps samples with fraction of Ns above specified threshold.")

    args.add_argument("--reference", type=str, help="If path to reference specified (FASTA), then whole genome will be written.")

    group = args.add_mutually_exclusive_group()

    group.add_argument("--include")
    group.add_argument("--exclude")

    args.add_argument("--with-stats", help="If a path is specified, then position of the outputed SNPs is stored in this file. ")

    return args.parse_args()

def main():
    """
    Process VCF files and merge them into a single fasta file.
    """

    logging.basicConfig(level=logging.INFO)

    args = get_args()
    contigs = list()

    sample_stats = dict()

    # All positions available for analysis.
    avail_pos = dict()
    # Stats about each position in each chromosome.
    pos_stats = dict()
    # Cached version of the data.
    vcf_data = dict()
    mixtures = dict()

    empty_tree = FastRBTree()

    exclude = False
    include = False

    if args.reference:
        ref_seq = OrderedDict()
        with open(args.reference) as fp:
            for record in SeqIO.parse(fp, "fasta"):
                ref_seq[record.id] = str(record.seq)

        args.reference = ref_seq

    if args.exclude or args.include:
        pos = {}
        chr_pos = []
        bed_file = args.include if args.include is not None else args.exclude

        with open(bed_file) as fp:
            for line in fp:
                data = line.strip().split("\t")

                chr_pos += [ (i, False,) for i in xrange(int(data[1]), int(data[2]) + 1)]

                if data[0] not in pos:
                    pos[data[0]] = []

                pos[data[0]] += chr_pos


        pos = {chrom: FastRBTree(l) for chrom, l in pos.items()}

        if args.include:
            include = pos
        else:
            exclude = pos


    if args.directory is not None and args.input is None:
        args.input = glob.glob(os.path.join(args.directory, "*.vcf"))

    # First pass to get the references and the positions to be analysed.
    for vcf_in in args.input:
        sample_name, _ = os.path.splitext(os.path.basename(vcf_in))
        vcf_data[vcf_in] = list()
        reader = vcf.Reader(filename=vcf_in)

        for record in reader:
            if include and include.get(record.CHROM, empty_tree).get(record.POS, True) or exclude and not exclude.get(record.CHROM, empty_tree).get(record.POS, True):
                continue

            vcf_data[vcf_in].append(record)

            if record.CHROM not in contigs:
                contigs.append(record.CHROM)
                avail_pos[record.CHROM] = FastRBTree()
                mixtures[record.CHROM] = {}
                sample_stats[record.CHROM] = {}

            if sample_name not in mixtures[record.CHROM]:
                mixtures[record.CHROM][sample_name] = FastRBTree()

            if sample_name not in sample_stats[record.CHROM]:
                sample_stats[record.CHROM][sample_name] = {}

            if record.FILTER == "PASS" or not record.FILTER:
                if record.is_snp:
                    if record.POS in avail_pos[record.CHROM] and avail_pos[record.CHROM][record.POS] != record.REF:
                        logging.critical("SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s", record.POS)
                        return 2

                    if record.CHROM not in pos_stats:
                        pos_stats[record.CHROM] = {}

                    avail_pos[record.CHROM].insert(record.POS, str(record.REF))
                    pos_stats[record.CHROM][record.POS] = {"N":0, "-": 0, "mut": 0}

            elif args.with_mixtures and record.is_snp:
                mix = get_mixture(record, args.with_mixtures)

                for ratio, code in mix.items():
                    for c in code:
                        avail_pos[record.CHROM].insert(record.POS, str(record.REF))
                        if record.CHROM not in pos_stats:
                            pos_stats[record.CHROM] = {}
                        pos_stats[record.CHROM][record.POS] = {"N": 0, "-": 0, "mut": 0}

                        if sample_name not in mixtures[record.CHROM]:
                            mixtures[record.CHROM][sample_name] = FastRBTree()

                        mixtures[record.CHROM][sample_name].insert(record.POS, c)


    all_data = { contig: {} for contig in contigs}
    samples = []

    for vcf_in in args.input:

        sample_seq = ""
        sample_name, _ = os.path.splitext(os.path.basename(vcf_in))
        samples.append(sample_name)

        # Initialise the data for this sample to be REF positions.
        for contig in contigs:
            all_data[contig][sample_name] = { pos: avail_pos[contig][pos] for pos in avail_pos[contig] }

#         reader = vcf.Reader(filename=vcf_in)
        for record in vcf_data[vcf_in]:
            # Array of filters that have been applied.
            filters = []

            # If position is our available position.
            if avail_pos.get(record.CHROM, empty_tree).get(record.POS, False):
                if record.FILTER == "PASS" or not record.FILTER:
                    if record.is_snp:
                        if len(record.ALT) > 1:
                            logging.info("POS %s passed filters but has multiple alleles. Inserting N")
                            all_data[record.CHROM][sample_name][record.POS] = "N"
                        else:
                            all_data[record.CHROM][sample_name][record.POS] = record.ALT[0].sequence
                            pos_stats[record.CHROM][record.POS]["mut"] += 1
                else:

                    # Currently we are only using first filter to call consensus.
                    extended_code = mixtures[record.CHROM][sample_name].get(record.POS, "N")

#                     extended_code = PHEFilterBase.call_concensus(record)

                    # Calculate the stats
                    if extended_code == "N":
                        pos_stats[record.CHROM][record.POS]["N"] += 1
                        pos_stats[record.CHROM][record.POS]["mut"] += 1

                        if "n_pos" not in sample_stats[record.CHROM][sample_name]:
                            sample_stats[record.CHROM][sample_name]["n_pos"] = []
                        sample_stats[record.CHROM][sample_name]["n_pos"].append(record.POS)

                    elif extended_code == "-":
                        pos_stats[record.CHROM][record.POS]["-"] += 1
#                     else:
#                         print "Good mixture %s: %i (%s)" % (sample_name, record.POS, extended_code)

                    # Save the extended code of the SNP.
                    all_data[record.CHROM][sample_name][record.POS] = extended_code
        del vcf_data[vcf_in]

    # Output the data to the fasta file.
    # The data is already aligned so simply output it.
    discarded = 0

    if args.reference:
        # These should be in the same order as the order in reference.
        contigs = args.reference.keys()

    if args.sample_Ns:
        delete_samples = []
        for contig in contigs:
            for sample in samples:

                # Skip if the contig not in sample_stats
                if contig not in sample_stats:
                    continue

                sample_n_ratio = float(len(sample_stats[contig][sample]["n_pos"])) / len(avail_pos[contig])
                if sample_n_ratio > args.sample_Ns:
                    for pos in sample_stats[contig][sample]["n_pos"]:
                        pos_stats[contig][pos]["N"] -= 1

                    logging.info("Removing %s due to high Ns in sample: %s", sample , sample_n_ratio)

                    delete_samples.append(sample)

        samples = [sample for sample in samples if sample not in delete_samples]
    snp_positions = []
    with open(args.out, "w") as fp:

        for sample in samples:
            sample_seq = ""
            for contig in contigs:
                if contig in avail_pos:
                    if args.reference:
                        positions = xrange(1, len(args.reference[contig]) + 1)
                    else:
                        positions = avail_pos[contig].keys()
                    for pos in positions:
                        if pos in avail_pos[contig]:
                            if not args.column_Ns or float(pos_stats[contig][pos]["N"]) / len(samples) < args.column_Ns and \
                                float(pos_stats[contig][pos]["-"]) / len(samples) < args.column_Ns:
                                sample_seq += all_data[contig][sample][pos]
                            else:
                                discarded += 1
                        elif args.reference:
                            sample_seq += args.reference[contig][pos - 1]
                elif args.reference:
                    sample_seq += args.reference[contig]



            fp.write(">%s\n%s\n" % (sample, sample_seq))
        # Do the same for reference data.
        ref_snps = ""

        for contig in contigs:
            if contig in avail_pos:
                if args.reference:
                    positions = xrange(1, len(args.reference[contig]) + 1)
                else:
                    positions = avail_pos[contig].keys()
                for pos in positions:
                    if pos in avail_pos[contig]:
                        if not args.column_Ns or float(pos_stats[contig][pos]["N"]) / len(samples) < args.column_Ns and \
                                float(pos_stats[contig][pos]["-"]) / len(samples) < args.column_Ns:

                            ref_snps += str(avail_pos[contig][pos])
                            snp_positions.append((contig, pos,))
                    elif args.reference:
                        ref_snps += args.reference[contig][pos - 1]
            elif args.reference:
                    ref_snps += args.reference[contig]

        fp.write(">reference\n%s\n" % ref_snps)

    if args.with_stats:
        with open(args.with_stats, "wb") as fp:
            for values in snp_positions:
                fp.write("%s\t%s\n" % (values[0], values[1]))
        plot_stats(pos_stats)
    # print_stats(sample_stats, pos_stats, total_vars=len(avail_pos[contig]))


    logging.info("Discarded total of %i poor quality columns", float(discarded) / len(args.input))
    return 0

if __name__ == '__main__':
    import time

#     with PyCallGraph(output=graphviz):
#     T0 = time.time()
    r = main()
#     T1 = time.time()

#     print "Time taken: %i" % (T1 - T0)
    exit(r)
