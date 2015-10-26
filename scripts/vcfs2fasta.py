'''
Merge SNP data from multiple VCF files into a single fasta file.

Created on 5 Oct 2015

@author: alex
'''
import argparse
import glob
import itertools
import os

from bintrees import FastRBTree
import vcf

from phe.variant_filters import IUPAC_CODES


def get_mixture(record, threshold):
    mixtures = {}
    try:
        if len(record.samples[0].data.AD) > 1:

            total_depth = sum(record.samples[0].data.AD)
            # Go over all combinations of touples.
            for comb in itertools.combinations(range(0, len(record.samples[0].data.AD)), 2):
                i = comb[0]
                j = comb[1]

                alleles = set()

                if 0 in comb:
                    alleles.add(str(record.REF))

                if i != 0:
                    alleles.add(str(record.ALT[i - 1]))
                    mixture = record.samples[0].data.AD[i]
                if j != 0:
                    alleles.add(str(record.ALT[j - 1]))
                    mixture = record.samples[0].data.AD[j]

                ratio = float(mixture) / total_depth
                if ratio == 1.0:
                    print "This is only designed for mixtures! %s %s %s" % (record, ratio, record.samples[0].data.AD)
                elif ratio >= threshold:
                    code = IUPAC_CODES[frozenset(alleles)]
                    if ratio not in mixtures:
                        mixtures[ratio] = []
                    mixtures[ratio].append(code)
    except AttributeError:
        mixtures = {}

    return mixtures

def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("--directory", "-d", help="Path to the directory with .vcf files.")
    args.add_argument("--input", "-i", type=str, nargs='+', help="List of VCF files to process.")
    args.add_argument("--out", "-o", required=True, help="Path to the output FASTA file.")

    args.add_argument("--with-mixtures", type=float)

    args.add_argument("--Ns", type=float)

    args.add_argument("--sample-Ns", type=float)

    return args.parse_args()

def main():
    """
    Process VCF files and merge them into a single fasta file.
    """

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

    if args.directory is not None and args.input is None:
        args.input = glob.glob(os.path.join(args.directory, "*.vcf"))

    # First pass to get the references and the positions to be analysed.
    for vcf_in in args.input:
        sample_name, _ = os.path.splitext(os.path.basename(vcf_in))
        vcf_data[vcf_in] = list()
        reader = vcf.Reader(filename=vcf_in)

        for record in reader:
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
                        print "SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s" % record.POS
                        return 2

                    if record.CHROM not in pos_stats:
                        pos_stats[record.CHROM] = {}

                    avail_pos[record.CHROM].insert(record.POS, str(record.REF))
                    pos_stats[record.CHROM][record.POS] = {"N":0, "-": 0}

            elif args.with_mixtures:
                mix = get_mixture(record, args.with_mixtures)

                for ratio, code in mix.items():
                    for c in code:
                        avail_pos[record.CHROM].insert(record.POS, str(record.REF))
                        if record.CHROM not in pos_stats:
                            pos_stats[record.CHROM] = {}
                        pos_stats[record.CHROM][record.POS] = {"N": 0, "-": 0}

                        if sample_name not in mixtures[record.CHROM]:
                            mixtures[record.CHROM][sample_name] = FastRBTree()

                        mixtures[record.CHROM][sample_name].insert(record.POS, c)


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

                    # Currently we are only using first filter to call consensus.
                    extended_code = mixtures[record.CHROM][sample_name].get(record.POS, "N")

#                     extended_code = PHEFilterBase.call_concensus(record)

                    # Calculate the stats
                    if extended_code == "N":
                        pos_stats[record.CHROM][record.POS]["N"] += 1

                        if "n_pos" not in sample_stats[record.CHROM][sample_name]:
                            sample_stats[record.CHROM][sample_name]["n_pos"] = []
                        sample_stats[record.CHROM][sample_name]["n_pos"].append(record.POS)

                    elif extended_code == "-":
                        pos_stats[record.CHROM][record.POS]["-"] += 1

                    # Save the extended code of the SNP.
                    all_data[record.CHROM][sample_name][record.POS] = extended_code

    # Output the data to the fasta file.
    # The data is already aligned so simply output it.
    discarded = 0

    if args.sample_Ns:
        delete_samples = []
        for contig in contigs:
            for sample in samples:
                sample_n_ratio = float(len(sample_stats[contig][sample]["n_pos"])) / len(avail_pos[contig])
                if sample_n_ratio > args.sample_Ns:
                    for pos in sample_stats[contig][sample]["n_pos"]:
                        pos_stats[contig][pos]["N"] -= 1

                    print "Removing %s due to high Ns in sample: %s" % (sample , sample_n_ratio)
                    delete_samples.append(sample)

        samples = [sample for sample in samples if sample not in delete_samples]

    with open(args.out, "w") as fp:

        for sample in samples:
            sample_seq = ""
            for contig in contigs:
                for pos in avail_pos[contig]:
                    if not args.Ns or float(pos_stats[contig][pos]["N"]) / len(samples) < args.Ns and \
                        float(pos_stats[contig][pos]["-"]) / len(samples) < args.Ns:
                        sample_seq += all_data[contig][sample][pos]
                    else:
                        discarded += 1


            fp.write(">%s\n%s\n" % (sample, sample_seq))
        # Do the same for reference data.
        ref_snps = ""
        for contig in contigs:
            for pos in avail_pos[contig]:
                if not args.Ns or float(pos_stats[contig][pos]["N"]) / len(samples) < args.Ns and \
                        float(pos_stats[contig][pos]["-"]) / len(samples) < args.Ns:
                    ref_snps += str(avail_pos[contig][pos])
        fp.write(">reference\n%s\n" % ref_snps)

    print("Discarded total of %i for poor quality columns" % (float(discarded) / len(args.input)))
    return 0

if __name__ == '__main__':
    import time

#     with PyCallGraph(output=graphviz):
#     T0 = time.time()
    r = main()
#     T1 = time.time()

#     print "Time taken: %i" % (T1 - T0)
    exit(r)
