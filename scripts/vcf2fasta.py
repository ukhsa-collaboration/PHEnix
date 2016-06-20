'''
Merge SNP data from multiple VCF files into a single fasta file.

:Date: 5 October, 2015
:Author: Alex Jironkin
'''

import argparse
from collections import OrderedDict
import glob
import itertools
import logging
import os
import shutil
import tempfile

from Bio import SeqIO
from bintrees import FastRBTree

from phe.utils.reader import ParallelVCFReader
from phe.variant_filters import IUPAC_CODES


# from phe.variant_filters import IUPAC_CODES
# Try importing the matplotlib and numpy for stats.
try:
    from matplotlib import pyplot as plt
    import numpy
    CAN_STATS = True
except ImportError:
    CAN_STATS = False

class BaseStats(object):
    """Simple class to keep statistics about a base."""
    def __init__(self, records=None):
        self.N = 0
        self.mut = 0
        self.gap = 0
        self.mix = 0
        self.total = 0
        self.NA = 0

    def __str__(self):
        return "N: %i, mut: %i, mix: %i, gap: %i, total: %i" % (self.N,
                                                                self.mut,
                                                                self.mix,
                                                                self.gap,
                                                                self.total)
    def __add__(self, other):
        self.mut += other.mut
        self.N += other.N
        self.gap += other.gap
        self.mix += other.mix
        self.total += other.total
        self.NA += other.NA

        return self

    def update(self, position_data, sample, reference):
        for k, v in position_data.iteritems():
            if k != sample:
                continue
            if v == "-":
                self.gap += 1
            elif v == "N":
                self.N += 1
            elif v in ["A", "C", "G", "T"]:
                self.mut += 1
            else:
                self.mix += 1
        self.total += 1

def _make_ref_insert(start, stop, reference, exclude):
    '''Create reference insert taking account exclude positions.'''
    if stop is None:
        return [c for i, c in enumerate(reference[start:]) if i + start not in exclude]
    else:
        return [c for i, c in enumerate(reference[start:stop - 1]) if i + start not in exclude]

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

def validate_record(record):
    """Validate the record."""
    if record.is_indel and not (record.is_uncallable or record.is_monomorphic) or len(record.REF) > 1:
        return False
    else:
        return True

def pick_best_records(records):
    """ Pick single record from multiple records for the same position.

    Parameters
    ----------
    records: dict
        Dictionary of lists containing records for the samples (1-many).

    Returns
    -------
    dict: Dictionary with a 1-1 mapping of samples to records.
    """

    final_selection = {}
    for k, v in records.iteritems():
        if len(v) == 1:
            if validate_record(v[0]):
                final_selection[k] = v[0]
            else:
                logging.debug("Discarding %s:%s from %s", v[0].CHROM, v[0].POS, k)
        else:
            logging.debug("Resolving multi-record ambiguity for %s (%s)", k, ",".join(str(r) for r in v))
            r = None
            for record in v:
                if not validate_record(record):
                    continue
                elif r is None:
                    r = record
                else:
                    if len(record.FILTER) > len(r.FILTER):
                        # Pick only those that have more failed filters, than current.
                        r = record
            if r is None:
                logging.error("Should have picked something, but everything failed.")
            final_selection[k] = r

    return final_selection

def get_mixture(record, threshold):
    """Generate proper UIPAC letter for mixture.

    Mixture is determined for any positions where ALT depth/total depth
    is **>=** than ``threshold``.

    Parameters
    ----------
    record: :py:class:`vcf._Record`
        Record to assign mixtures to.
    threshold: float
        Threshold above which mixtures are calculated
    """

    mixture = "N"
    try:
        if len(record.samples[0].data.AD) > 1:

            total_depth = sum(record.samples[0].data.AD)
            # Go over all combinations of touples.
            alleles = list()
            for i in range(0, len(record.samples[0].data.AD)):

                mixture = record.samples[0].data.AD[i]

                ratio = float(mixture) / total_depth
                if ratio == 1.0:
                    logging.debug("This is only designed for mixtures! %s %s %s %s", record, ratio, record.samples[0].data.AD, record.FILTER)

                elif ratio >= threshold:
                    if i == 0:
                        alleles.append(str(record.REF))
                    else:
                        alleles.append(str(record.ALT[i - 1]))
                else:
                    logging.debug("Discarding allele %s to be below threshold %s", record.ALT[i - 1], ratio)
            try:
                mixture = IUPAC_CODES[frozenset(alleles)]

            except KeyError:
                # This happens when we have empty set of alleles.
                #    So we don't know what to put in there -> N.
                logging.debug("Could not retrieve IUPAC code for %s from %s", alleles, record)
                mixture = "N"

    except AttributeError:
        mixture = "N"

    return mixture

def is_above_min_depth(record):

    above = True
    for f in record.FILTER:
        if "min_depth" in f:
            above = False
            break
    return above

def print_stats(stats, pos_stats, total_vars):
    for contig in stats:
        for sample, info in stats[contig].items():
            print "%s,%i,%i" % (sample, len(info.get("n_pos", [])), total_vars)

    for contig in stats:
        for pos, info in pos_stats[contig].iteritems():
            print "%s,%i,%i,%i,%i" % (contig, pos, info.get("N", "NA"), info.get("-", "NA"), info.get("mut", "NA"))

def get_desc():
    return "Combine multiple VCFs into a single FASTA file."

def get_args():

    def positive_float(value):
        """Make type for float 0<x<1."""
        x = float(value)
        if not 0.0 <= x <= 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % x)
        return x

    args = argparse.ArgumentParser(description=get_desc())

    group = args.add_mutually_exclusive_group(required=True)
    group.add_argument("--directory", "-d", help="Path to the directory with .vcf files.")
    group.add_argument("--input", "-i", type=str, nargs='+', help="List of VCF files to process.")

    args.add_argument("--regexp", type=str, help="Regular expression for finding VCFs in a directory.")

    args.add_argument("--out", "-o", required=True, help="Path to the output FASTA file.")

    args.add_argument("--with-mixtures", type=positive_float, help="Specify this option with a threshold to output mixtures above this threshold.")

    args.add_argument("--column-Ns", type=positive_float, help="Keeps columns with fraction of Ns below specified threshold.")
    args.add_argument("--column-gaps", type=positive_float, help="Keeps columns with fraction of Ns below specified threshold.")

    args.add_argument("--sample-Ns", type=positive_float, help="Keeps samples with fraction of Ns below specified threshold.")
    args.add_argument("--sample-gaps", type=positive_float, help="Keeps samples with fraction of gaps below specified threshold.")

    args.add_argument("--reference", type=str, help="If path to reference specified (FASTA), then whole genome will be written.")

    group = args.add_mutually_exclusive_group()

    group.add_argument("--include", help="Only include positions in BED file in the FASTA")
    group.add_argument("--exclude", help="Exclude any positions specified in the BED file.")

    args.add_argument("--with-stats", help="If a path is specified, then position of the outputed SNPs is stored in this file. Requires mumpy and matplotlib.")
#     args.add_argument("--plots-dir", default="plots", help="Where to write summary plots on SNPs extracted. Requires mumpy and matplotlib.")

    return args

def main(args):
    """
    Process VCF files and merge them into a single fasta file.
    """

    contigs = list()

    empty_tree = FastRBTree()

    exclude = {}
    include = {}
    out_dir = os.path.join(os.path.dirname(args["out"]), "tmp")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if args["reference"]:
        ref_seq = OrderedDict()
        with open(args["reference"]) as fp:
            for record in SeqIO.parse(fp, "fasta"):
                ref_seq[record.id] = list(record.seq)

        args["reference"] = ref_seq

    if args["exclude"] or args["include"]:
        pos = {}
        chr_pos = []
        bed_file = args["include"] if args["include"] is not None else args["exclude"]

        with open(bed_file) as fp:
            for line in fp:
                data = line.strip().split("\t")

                chr_pos += [ (i, False,) for i in xrange(int(data[1]), int(data[2]) + 1)]

                if data[0] not in pos:
                    pos[data[0]] = []

                pos[data[0]] += chr_pos

        pos = {chrom: FastRBTree(l) for chrom, l in pos.items()}

        if args["include"]:
            include = pos
        else:
            exclude = pos


    if args["directory"] is not None and args["input"] is None:
        regexp = args["regexp"] if args["regexp"] else "*.vcf"
        args["input"] = glob.glob(os.path.join(args["directory"], regexp))

    if not args["input"]:
        logging.warn("No VCFs found.")
        return 0


    # If we can stats and asked to stats, then output the data
    if args["with_stats"] is not None:
        args["with_stats"] = open(args["with_stats"], "wb")
        args["with_stats"].write("contig,position,mutations,n_frac,n_gaps\n")


    parallel_reader = ParallelVCFReader(args["input"])

    sample_seqs = { sample_name: tempfile.NamedTemporaryFile(prefix=sample_name, dir=out_dir) for sample_name in parallel_reader.get_samples() }
    sample_seqs["reference"] = tempfile.NamedTemporaryFile(prefix="reference", dir=out_dir)

    samples = parallel_reader.get_samples() + ["reference"]
    sample_stats = {sample: BaseStats() for sample in samples }
    last_base = 0

    for chrom, pos, records in parallel_reader:

        final_records = pick_best_records(records)
        reference = [ record.REF for record in final_records.itervalues() if record.REF != "N"]
        valid = not reference or reference.count(reference[0]) == len(reference)

        # Make sure reference is the same across all samples.
        assert valid, "Position %s is not valid as multiple references found: %s" % (pos, reference)

        if not reference:
            continue
        else:
            reference = reference[0]

        # SKIP (or include) any pre-specified regions.
        if include and pos not in include.get(chrom, empty_tree) or exclude and pos in exclude.get(chrom, empty_tree):
            continue

        position_data = {"reference": str(reference), "stats": BaseStats()}

        for sample_name, record in final_records.iteritems():

            position_data["stats"].total += 1

            # IF this is uncallable genotype, add gap "-"
            if record.is_uncallable:
                # TODO: Mentioned in issue: #7(gitlab)
                position_data[sample_name] = "-"

                # Update stats
                position_data["stats"].gap += 1


            elif not record.FILTER:
                # If filter PASSED!
                # Make sure the reference base is the same. Maybe a vcf from different species snuck in here?!
                assert str(record.REF) == position_data["reference"] or str(record.REF) == 'N' or position_data["reference"] == 'N', "SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s in %s (%s, %s)" % (record.POS, sample_name, str(record.REF), position_data["reference"])
                # update position_data['reference'] to a real base if possible
                if position_data['reference'] == 'N' and str(record.REF) != 'N':
                    position_data['reference'] = str(record.REF)
                if record.is_snp:
                    if len(record.ALT) > 1:
                        logging.info("POS %s passed filters but has multiple alleles REF: %s, ALT: %s. Inserting N", record.POS, str(record.REF), str(record.ALT))
                        position_data[sample_name] = "N"
                        position_data["stats"].N += 1

                    else:
                        position_data[sample_name] = str(record.ALT[0])

                        position_data["stats"].mut += 1

            # Filter(s) failed
            elif record.is_snp and is_above_min_depth(record):
                if args["with_mixtures"]:
                    extended_code = get_mixture(record, args["with_mixtures"])
                else:
                    extended_code = "N"

                if extended_code == "N":
                    position_data["stats"].N += 1
                elif extended_code in ["A", "C", "G", "T"]:
                    position_data["stats"].mut += 1
                else:
                    position_data["stats"].mix += 1

                position_data[sample_name] = extended_code

            else:
                # filter fail; code as N for consistency
                position_data[sample_name] = "N"
                position_data["stats"].N += 1

            # Filter columns when threashold reaches user specified value.
            if isinstance(args["column_Ns"], float) and float(position_data["stats"].N) / len(args["input"]) > args["column_Ns"]:
                break
#                 del position_data[sample_name]

            if isinstance(args["column_gaps"], float) and float(position_data["stats"].gap) / len(args["input"]) > args["column_gaps"]:
                break
#                 del position_data[sample_name]

        else:
            if args["reference"]:
                seq = _make_ref_insert(last_base, pos, args["reference"][chrom], exclude.get(chrom, empty_tree))
                for sample in samples:
#                     sample_seqs[sample] += seq
                    sample_seqs[sample].write(''.join(seq))

            for i, sample_name in enumerate(samples):
                sample_base = position_data.get(sample_name, reference)

#                 sample_seqs[sample_name] += [sample_base]
                sample_seqs[sample_name].write(sample_base)
                sample_stats[sample_name].update(position_data, sample_name, reference)

            if args["with_stats"] is not None:
                args["with_stats"].write("%s,%i,%0.5f,%0.5f,%0.5f\n" % (chrom,
                                             pos,
                                             float(position_data["stats"].mut) / len(args["input"]),
                                             float(position_data["stats"].N) / len(args["input"]),
                                             float(position_data["stats"].gap) / len(args["input"]))
                         )

            last_base = pos

    # Fill from last snp to the end of reference.
    # FIXME: A little naughty to use chrom outside the loop!
    if args["reference"]:
        seq = _make_ref_insert(last_base, None, args["reference"][chrom], exclude.get(chrom, empty_tree))
        for sample in samples:
#             sample_seqs[sample] += seq
            sample_seqs[sample].write(''.join(seq))

    sample_seqs["reference"].seek(0)
    reference = sample_seqs["reference"].next()
    sample_seqs["reference"].close()
    del sample_seqs["reference"]

    # Exclude any samples with high Ns or gaps
    if isinstance(args["sample_Ns"], float):
        for sample_name in samples:
            if sample_name == "reference":
                continue
            n_fraction = float(sample_stats[sample_name].N) / sample_stats[sample_name].total
            if n_fraction > args["sample_Ns"]:
                logging.info("Removing %s due to high sample Ns fraction %s", sample_name, n_fraction)

                sample_seqs[sample_name].close()
                del sample_seqs[sample_name]
                del sample_stats[sample_name]

    # Exclude any samples with high gap fraction.
    if isinstance(args["sample_gaps"], float):
        for sample_name in samples:
            if sample_name == "reference":
                continue

            gap_fractoin = float(sample_stats[sample_name].gap) / sample_stats[sample_name].total
            if gap_fractoin > args["sample_gaps"]:
                logging.info("Removing %s due to high sample gaps fraction %s", sample_name, gap_fractoin)

                sample_seqs[sample_name].close()
                del sample_seqs[sample_name]
                del sample_stats[sample_name]

    try:
        assert len(sample_seqs) > 0, "All samples have been filtered out."

        with open(args["out"], "w") as fp:
            fp.write(">reference\n%s\n" % reference)
            reference_length = len(reference)
            del reference

            for sample_name, tmp_iter in sample_seqs.iteritems():
                tmp_iter.seek(0)
                # These are dumped as single long string of data. Calling next() should read it all.
                snp_sequence = tmp_iter.next()
                assert len(snp_sequence) == reference_length, "Sample %s has length %s, but should be %s (reference)" % (sample_name, len(snp_sequence), reference_length)

                fp.write(">%s\n%s\n" % (sample_name, ''.join(snp_sequence)))
    except AssertionError as e:
        logging.error(e.message)

        # Need to delete the malformed file.
        os.unlink(args["out"])

    finally:
        # Close all the tmp handles.
        for tmp_iter in sample_seqs.itervalues():
            tmp_iter.close()
        shutil.rmtree(out_dir)

        if args["with_stats"] is not None:
            args["with_stats"].close()

    # Compute the stats.
    for sample in sample_stats:
        if sample != "reference":
            print "%s\t%s" % (sample, str(sample_stats[sample]))

#         if CAN_STATS:
#             plot_stats(avail_pos, len(samples) - 1, plots_dir=os.path.abspath(args["plots_dir"]))

    return 0

if __name__ == '__main__':
    exit(main(vars(get_args().parse_args())))
