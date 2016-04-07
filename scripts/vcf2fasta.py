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
import vcf

from phe.variant_filters import IUPAC_CODES


# Try importing the matplotlib and numpy for stats.
try:
    from matplotlib import pyplot as plt
    import numpy
    can_stats = True
except ImportError:
    can_stats = False

def is_uncallable(record):

    uncall = False
    try:
        if record.samples[0].data.GT in ("./.", None):
            uncall = True
    except:
        uncall = None

    if record.FILTER is not None and "LowQual" in record.FILTER:
        uncall = True

    return uncall


class base_stats(object):
    def __init__(self):
        self.N = 0
        self.mut = 0
        self.gap = 0
        self.mix = 0
        self.total = 0

    def __str__(self):
        return "N: %i, mut: %i, mix: %i, gap: %i, total: %i" % (self.N, self.mut, self.mix, self.gap, self.total)

def get_sample_stats(all_positions, samples):
    sample_stats = {sample: base_stats() for sample in samples }
    for positions in all_positions.itervalues():
        for position in positions:
            for sample in samples:
                base = positions[position].get(sample)

                if base == "-":
                    sample_stats[sample].gap += 1
                elif base == "N":
                    sample_stats[sample].N += 1
                elif base is not None and base != positions[position].get("reference"):
                    sample_stats[sample].mut += 1

    return sample_stats

def calculate_dist_matrix():
    return None

def plot_stats(pos_stats, total_samples, plots_dir="plots", discarded={}):
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    for contig in pos_stats:
        plt.style.use('ggplot')

        x = numpy.array([pos for pos in pos_stats[contig] if pos not in discarded.get(contig, [])])
        y = numpy.array([ float(pos_stats[contig][pos]["stats"].mut) / total_samples for pos in pos_stats[contig] if pos not in discarded.get(contig, []) ])

        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
        f.set_size_inches(12, 15)
        ax1.plot(x, y, 'ro')
        ax1.set_title("Fraction of samples with SNPs")
        plt.ylim(0, 1.1)

        y = numpy.array([ float(pos_stats[contig][pos]["stats"].N) / total_samples for pos in pos_stats[contig] if pos not in discarded.get(contig, [])])
        ax2.plot(x, y, 'bo')
        ax2.set_title("Fraction of samples with Ns")

        y = numpy.array([ float(pos_stats[contig][pos]["stats"].mix) / total_samples for pos in pos_stats[contig] if pos not in discarded.get(contig, [])])
        ax3.plot(x, y, 'go')
        ax3.set_title("Fraction of samples with mixed bases")

        y = numpy.array([ float(pos_stats[contig][pos]["stats"].gap) / total_samples for pos in pos_stats[contig] if pos not in discarded.get(contig, [])])
        ax4.plot(x, y, 'yo')
        ax4.set_title("Fraction of samples with uncallable genotype (gap)")

        contig = contig.replace("/", "-")
        plt.savefig(os.path.join(plots_dir, "%s.png" % contig), dpi=100)

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

    args.add_argument("--regexp", type=str, help="Regular expression for finding VCFs in a directory.")

    args.add_argument("--out", "-o", required=True, help="Path to the output FASTA file.")

    args.add_argument("--with-mixtures", type=float, help="Specify this option with a threshold to output mixtures above this threshold.")

    args.add_argument("--column-Ns", type=float, help="Keeps columns with fraction of Ns above specified threshold.")
    args.add_argument("--column-gaps", type=float, help="Keeps columns with fraction of Ns above specified threshold.")

    args.add_argument("--sample-Ns", type=float, help="Keeps samples with fraction of Ns above specified threshold.")

    args.add_argument("--reference", type=str, help="If path to reference specified (FASTA), then whole genome will be written.")

    group = args.add_mutually_exclusive_group()

    group.add_argument("--include", help="Only include positions in BED file in the FASTA")
    group.add_argument("--exclude", help="Exclude any positions specified in the BED file.")

    args.add_argument("--with-stats", help="If a path is specified, then position of the outputed SNPs is stored in this file. Requires mumpy and matplotlib.")
    args.add_argument("--plots-dir", default="plots", help="Where to write summary plots on SNPs extracted. Requires mumpy and matplotlib.")

    args.add_argument("--with-dist-mat", help="If this option is specified, then distance matrix is calculated for the samples.")
    args.add_argument("--count-dist-gaps", action="store_true", help="Counr gaps as valid character in distance matrix.")
    args.add_argument("--count-dist-Ns", action="store_true", help="Counr Ns as valid character in distance matrix.")

    args.add_argument("--debug", action="store_true", help="More verbose logging (default: turned off).")

    return args

def main():
    """
    Process VCF files and merge them into a single fasta file.
    """

    args = get_args().parse_args()

    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO)

    contigs = list()

    samples = list()
    valid_chars = ["A", "C", "G", "T"]

    if args.count_dist_gaps:
        valid_chars.append("-")

    if args.count_dist_Ns:
        valid_chars.append("N")

    # All positions available for analysis.
    avail_pos = dict()

    empty_tree = FastRBTree()

    exclude = {}
    include = {}

    if args.reference:
        ref_seq = OrderedDict()
        with open(args.reference) as fp:
            for record in SeqIO.parse(fp, "fasta"):
                ref_seq[record.id] = list(record.seq)

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
        regexp = args.regexp if args.regexp else "*.vcf"
        args.input = glob.glob(os.path.join(args.directory, regexp))
    else:
        args.input = [args.input]

    if not args.input:
        logging.warn("No VCFs found.")
        return 0

    # First pass to get the references and the positions to be analysed.
    for vcf_in in args.input:
        sample_name, _ = os.path.splitext(os.path.basename(vcf_in))
        samples.append(sample_name)

        reader = vcf.Reader(filename=vcf_in)

        # Go over every position in the reader.
        for record in reader:

            # Inject a property about uncallable genotypes into record.
            record.__setattr__("is_uncallable", is_uncallable(record))  # is_uncallable = types.MethodType(is_uncallable, record)

            # SKIP indels, if not handled then can cause REF base to be >1
            if record.is_indel and not (record.is_uncallable or record.is_monomorphic) or len(record.REF) > 1:
            # if len(record.REF) > 1:
                # print "%s\t%s\t%s\t%s\t%s" % (sample_name,record.POS,-1,record.FILTER,record)
                continue
                # if record.is_deletion and not record.is_uncallable:
                #    continue

            # SKIP (or include) any pre-specified regions.
            if include and record.POS not in include.get(record.CHROM, empty_tree) or exclude and record.POS in exclude.get(record.CHROM, empty_tree):
#             if include.get(record.CHROM, empty_tree).get(record.POS, False) or \
#                 not exclude.get(record.CHROM, empty_tree).get(record.POS, True):
                continue

            # Setup the RB tree for contigs not in the data structure.
            if record.CHROM not in contigs:
                contigs.append(record.CHROM)
                avail_pos[record.CHROM] = FastRBTree()

            # Setup the position data to contain reference and stats.
            if avail_pos[record.CHROM].get(record.POS, None) is None:
                avail_pos[record.CHROM].insert(record.POS, {"reference": str(record.REF),
                                                            "stats": base_stats()})

            position_data = avail_pos[record.CHROM].get(record.POS)

            assert len(position_data["reference"]) == 1, "Reference base must be singluar: in %s found %s @ %s" % (position_data["reference"], sample_name, record.POS)
            where = 0
            # IF this is uncallable genotype, add gap "-"
            if record.is_uncallable:
                # TODO: Mentioned in issue: #7(gitlab)
                position_data[sample_name] = "-"

                # Update stats
                position_data["stats"].N += 1
                where = 1

            elif not record.FILTER:
                # If filter PASSED!
                # Make sure the reference base is the same. Maybe a vcf from different species snuck in here?!
                assert str(record.REF) == position_data["reference"], "SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s in %s (%s, %s)" % (record.POS, vcf_in, str(record.REF), position_data["reference"])
                if record.is_snp:
                    if len(record.ALT) > 1:
                        logging.info("POS %s passed filters but has multiple alleles. Inserting N")
                        position_data[sample_name] = "N"
                        position_data["stats"].N += 1

                    else:
                        position_data[sample_name] = str(record.ALT[0])

                        position_data["stats"].mut += 1
                where = 2
            # Filter(s) failed
            elif record.is_snp:
                # mix = get_mixture(record, args.with_mixtures)
                # Currently we are only using first filter to call consensus.
                extended_code = "N"

                if extended_code == "N":
                    position_data["stats"].N += 1

                position_data[sample_name] = extended_code
                where = 3
            else:
                # filter fail; code as N for consistency
                position_data[sample_name] = "N"
                position_data["stats"].N += 1
                where = 4
            # print "%s\t%s\t%s\t%s\t%s" % (sample_name,record.POS,where,record.FILTER,record)
            # For reference we always want to use all data.
            if args.reference:
                continue

            # Filter columns when threashold reaches user specified value.
            if isinstance(args.column_Ns, float) and float(position_data["stats"].N) / len(args.input) > args.column_Ns:
                avail_pos[record.CHROM].remove(record.POS)

                # print "excluding %s" % record.POS
                if record.CHROM not in exclude:
                    exclude[record.CHROM] = FastRBTree()
                exclude[record.CHROM].insert(record.POS, False)

            if isinstance(args.column_gaps, float) and float(position_data["stats"].gap) / len(args.input) > args.column_gaps:
                avail_pos[record.CHROM].remove(record.POS)

                if record.CHROM not in exclude:
                    exclude[record.CHROM] = FastRBTree()
                exclude[record.CHROM].insert(record.POS, False)

    # Compute per sample statistics.
    sample_stats = get_sample_stats(avail_pos, samples)

    # Exclude any samples with high Ns or gaps
    if isinstance(args.sample_Ns, float):
        ss = []
        for sample_name in samples:
            total_positions = 0
            for contig in contigs:
                total_positions += len(avail_pos[contig])
            if sample_stats[sample_name].N / total_positions <= args.sample_Ns:
                ss.append(sample_name)
            else:
                pass
                # print "EXCLUDING: %s for high N fraction: %s" % (sample_name, sample_stats[sample_name].N / total_positions)

        samples = ss

    # ALWAYS APPEND reference
    samples.append("reference")
    dist_mat = {}
    sample_seqs = { sample_name: [] for sample_name in samples }
    c = 0

    if args.with_dist_mat:
        for i, sample_1 in enumerate(samples):
            dist_mat[sample_1] = {}
            for j, sample_2 in enumerate(samples):
                if j < i:
                    continue
                dist_mat[sample_1][sample_2] = 0

    # For each contig concatinate sequences.
    for contig in contigs:

        # if contig is not in the avail pos then concatinate the whole reference.
        if args.reference and contig not in avail_pos:
            for sample in samples:
                sample_seqs[sample] += args.reference[contig]
            continue

        last_base = 0
        for pos in avail_pos[contig]:
            c += 1
            if args.reference:
                seq = args.reference[contig][last_base:pos - 1]
                for sample in samples:
                    sample_seqs[sample] += seq

            bases = set()
            ref_base = avail_pos[contig][pos].get("reference")
            # Position has been seen or no reference available.
            for i, sample in enumerate(samples):

                sample_base = avail_pos[contig][pos].get(sample, ref_base)

                sample_seqs[sample] += [sample_base]
                bases.add(sample_base)

                # If we don't need distance matrix, then continue from top.
                if not args.with_dist_mat or sample_base.upper() not in valid_chars:
                    continue

                for j, sample_2 in enumerate(samples):
                    if j <= i:
                        continue

                    s2_base = avail_pos[contig][pos].get(sample_2, ref_base)

                    if sample_base != s2_base and s2_base.upper() in valid_chars:
                        dist_mat[sample][sample_2] += 1

            # Do the internal check that positions have at least 2 different characters.
            # assert len(bases) > 1, "Internal consustency check failed for position %s bases: %s" % (pos, bases)
            # removed the check because when removing a sample based on sample-Ns can lead to non SNP bases in column.

            # Keep track of the last processed position.
            last_base = pos

        # Fill from last snp to the end of reference.
        if args.reference:
            seq = args.reference[contig][last_base:]
            for sample in samples:
                sample_seqs[sample] += seq

    # Write the sequences out.
    with open(args.out, "w") as fp:
        for sample in sample_seqs:
            fp.write(">%s\n%s\n" % (sample, ''.join(sample_seqs[sample])))

    # Compute the stats.
    for sample in sample_stats:
        total_positions = 0
        for contig in contigs:
            total_positions += len(avail_pos[contig])
        sample_stats[sample].total = total_positions
        print "%s\t%s" % (sample, str(sample_stats[sample]))

    if args.with_dist_mat:
        with open(args.with_dist_mat, "wb") as fp:
            fp.write(",%s\n" % ",".join(samples))
            for i, sample_1 in enumerate(samples):
                row = "%s" % sample_1
                for j, sample_2 in enumerate(samples):
                    if j < i:
                        dist = dist_mat[sample_2][sample_1]
                    else:
                        dist = dist_mat[sample_1][sample_2]
                    row += ",%i" % dist
                fp.write("%s\n" % row)

    # If we can stats and asked to stats, then output the data
    if args.with_stats:
        with open(args.with_stats, "wb") as fp:
            fp.write("contig,position,mutations,n_frac,n_gaps\n")
            for contig in contigs:
                for pos in avail_pos[contig]:
                    position_data = avail_pos[contig][pos]
                    fp.write("%s,%i,%0.5f,%0.5f,%0.5f\n" % (contig,
                                                 pos,
                                                 float(position_data["stats"].mut) / len(args.input),
                                                 float(position_data["stats"].N) / len(args.input),
                                                 float(position_data["stats"].gap) / len(args.input))
                             )
        if can_stats:
            plot_stats(avail_pos, len(samples) - 1, plots_dir=os.path.abspath(args.plots_dir))

    return 0

if __name__ == '__main__':
    exit(main())
