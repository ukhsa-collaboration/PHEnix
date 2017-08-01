'''
Merge SNP data from multiple VCF files into a single fasta file.

:Date: 5 October, 2015
:Author: Alex Jironkin
'''

import sys
import argparse
from collections import OrderedDict
import glob
import logging
import os
import shutil
import tempfile
import vcf
from math import sqrt

from Bio import SeqIO

from phe.utils import is_uncallable
from phe.variant_filters import IUPAC_CODES

# --------------------------------------------------------------------------------------------------

def validate_record(record):
    """
    Validate the record.

    Parameters:
    -----------
    record: obj
        a record from the vcf library

    Returns:
    --------
    True or False
    """

    if record.is_indel and not (is_uncallable(record) or record.is_monomorphic) or len(record.REF) > 1:
        return False
    else:
        return True

# --------------------------------------------------------------------------------------------------

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
                    logging.debug("This is only designed for mixtures! %s %s %s %s",
                                  record,
                                  ratio,
                                  record.samples[0].data.AD,
                                  record.FILTER)

                elif ratio >= threshold:
                    if i == 0:
                        alleles.append(str(record.REF))
                    else:
                        alleles.append(str(record.ALT[i - 1]))
                else:
                    logging.debug("Discarding allele %s to be below threshold %s",
                                  record.ALT[i - 1],
                                  ratio)
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

# --------------------------------------------------------------------------------------------------

def is_above_min_depth(record):
    """
    Check if the record is above the min depth threshold.

    Parameters:
    -----------
    record: obj
        a record from the vcf library

    Returns:
    --------
    above: boolean
        True if record is above min depth
    """

    above = True
    for f in record.FILTER:
        if "min_depth" in f:
            above = False
            break
    return above

# --------------------------------------------------------------------------------------------------

def get_desc():
    """
    Return the description as a string

    Parameters:
    -----------

    Returns:
    --------
    "Combine multiple VCFs into a single FASTA file."

    """

    return "Combine multiple VCFs into a single FASTA file."

# --------------------------------------------------------------------------------------------------

def get_args():
    """
    Get the arguments.

    Parameters:
    -----------
    nil

    Returns:
    --------
    args: obj
        argument object
    """

    def positive_float(value):
        """Make type for float 0<x<1."""
        x = float(value)
        if not 0.0 <= x <= 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % x)
        return x

    args = argparse.ArgumentParser(description=get_desc())

    group_a = args.add_mutually_exclusive_group(required=True)
    group_a.add_argument("--directory",
                         "-d",
                         help="Path to the directory with .vcf files.")
    group_a.add_argument("--input",
                         "-i",
                         type=str,
                         nargs='+',
                         help="List of VCF files to process.")

    args.add_argument("--regexp",
                      type=str,
                      help="Regular expression for finding VCFs in a directory.")

    args.add_argument("--out",
                      "-o",
                      required=True,
                      help="Path to the output FASTA file.")

    args.add_argument("--with-mixtures",
                      type=positive_float,
                      help="Specify this option with a threshold to output mixtures above this threshold.")

    args.add_argument("--column-Ns",
                      type=positive_float,
                      help="Keeps columns with fraction of Ns below specified threshold.")
    args.add_argument("--column-gaps",
                      type=positive_float,
                      help="Keeps columns with fraction of Ns below specified threshold.")

    args.add_argument("--sample-Ns",
                      help="Keeps samples with fraction of Ns below specified threshold or put 'auto'." + \
                           "Fraction expressed as fraction of genome. Requires --reflength or --reference.")
    args.add_argument("--sample-gaps",
                      help="Keeps samples with fraction of gaps below specified threshold or put 'auto'." + \
                           "Fraction expressed as fraction of genome. Requires --reflength or --reference.")
    args.add_argument("--sample-Ns-gaps-auto-factor",
                      default=2.0,
                      type=float,
                      help="When using 'auto' option for --sample-gaps or --sample-Ns, remove sample that have" + \
                           "gaps or Ns this many times above the stddev of all samples. [Default: 2.0]")

    group_b = args.add_mutually_exclusive_group()
    group_b.add_argument("--reference",
                         type=str,
                         help="If path to reference specified (FASTA), then whole genome will be written to alignment.")
    group_b.add_argument("--remove-invariant-npos",
                         action='store_true',
                         help="Remove all positions that invariant apart from N positions.")

    args.add_argument("--reflength",
                      help="Length of reference. Either as int or can be worked out from fasta file. Ignored if --reference is used.")

    group_c = args.add_mutually_exclusive_group()
    group_c.add_argument("--include",
                         help="Only include positions in BED file in the FASTA")
    group_c.add_argument("--exclude",
                         help="Exclude any positions specified in the BED file.")

    args.add_argument("--with-stats",
                      help="If a path is specified, then position of the outputed SNPs is stored in this file.")

    return args

# --------------------------------------------------------------------------------------------------

def main(args):
    """
    Process VCF files and merge them into a single fasta file.

    Parameters:
    -----------
    args: obj
        argument object as dict

    Returns:
    --------
    0
    """

    # do some additional args checking
    if args['sample_Ns']:
        if args['reference'] == None and args['reflength'] == None:
            logging.error("--reflength or --reference is REQUIRED when using sample-Ns filter.")
            return 1
        if args['sample_Ns'] != 'auto':
            try:
                args['sample_Ns'] = float(args['sample_Ns'])
                if not 0.0 <= args['sample_Ns'] <= 1.0:
                    raise TypeError
            except TypeError:
                logging.error("Please put either 'auto' or a float [0.0, 1.0] in sample-Ns.")
                return 1

    if args['sample_gaps']:
        if args['reference'] == None and args['reflength'] == None:
            logging.error("--reflength or --reference is REQUIRED when using sample-gaps filter.")
            return 1
        if args['sample_gaps'] != 'auto':
            try:
                args['sample_gaps'] = float(args['sample_gaps'])
                if not 0.0 <= args['sample_gaps'] <= 1.0:
                    raise TypeError
            except TypeError:
                logging.error("Please put either 'auto' or a float [0.0, 1.0] in sample-gaps.")
                return 1

    if args["reference"]:
        if args['column_Ns']:
            logging.error("Reference option is given. This means that all columns" + \
                          " will be in the alignment, so it's mutually exclusive with" + \
                          " column-Ns option. Please chose either, but not both.")
            return 1

        ref_seq = OrderedDict()
        with open(args["reference"]) as fp:
            for record in SeqIO.parse(fp, "fasta"):
                ref_seq[record.id] = list(record.seq)

        args["reference"] = ref_seq
        args['reflength'] = sum([len(x) for x in ref_seq.values()])

    elif args['reflength']:
        try:
            args['reflength'] = int(args['reflength'])
        except ValueError:
            with open(args["reflength"]) as fp:
                args['reflength'] = sum([len(line.strip()) if not line.startswith('>') else 0 for line in fp])
    else:
        pass

    if args["directory"] is not None and args["input"] is None:
        regexp = args["regexp"] if args["regexp"] else "*.vcf"
        args["input"] = glob.glob(os.path.join(args["directory"], regexp))

    if not args["input"]:
        logging.warn("No VCFs found.")
        return 0

    nof_vcfs = len(args['input'])
    logging.info("%i vcf files found.", nof_vcfs)

    all_contig_data = {}

    cnt = 0
    # iterate over all vcfs files
    for vcf_file in args["input"]:
        cnt += 1

        reader = vcf.Reader(filename=vcf_file)
        sample_name = reader.samples[0]
        for record in reader:

            # SKIP indels, if not handled then can cause REF base to be >1
            if validate_record(record) == False:
                continue

            if all_contig_data.has_key(record.CHROM) == False:
                all_contig_data[record.CHROM] = {'reference': {}}

            data = all_contig_data[record.CHROM]

            if data.has_key(sample_name) == False:
                data[sample_name] = {}

            # IF this is uncallable genotype, add gap "-"
            if is_uncallable(record):
                try:
                    data[sample_name]['-'].add(record.POS)
                except KeyError:
                    data[sample_name]['-'] = set([record.POS])

            elif not record.FILTER:
                # If filter PASSED!
                if record.is_snp:
                    if len(record.ALT) > 1:
                        logging.info("POS %s passed filters but has multiple alleles REF: %s, ALT: %s. Inserting N",
                                     record.POS,
                                     str(record.REF),
                                     str(record.ALT))
                        try:
                            data[sample_name]['N'].add(record.POS)
                        except KeyError:
                            data[sample_name]['N'] = set([record.POS])
                    else:
                        try:
                            data[sample_name][str(record.ALT[0])].add(record.POS)
                        except KeyError:
                            data[sample_name][str(record.ALT[0])] = set([record.POS])

            # Filter(s) failed
            elif record.is_snp and is_above_min_depth(record):
                if args["with_mixtures"]:
                    extended_code = get_mixture(record, args["with_mixtures"])
                    if extended_code == str(record.REF):
                        # expressing this as a mixture would make it a ref base
                        # the right thing to do now is ignore the position completely
                        continue
                else:
                    extended_code = "N"
                try:
                    data[sample_name][extended_code].add(record.POS)
                except KeyError:
                    data[sample_name][extended_code] = set([record.POS])

            else:
                try:
                    data[sample_name]['N'].add(record.POS)
                except KeyError:
                    data[sample_name]['N'] = set([record.POS])

            # add all positions that are processed to the reference sets
            try:
                data['reference'][str(record.REF)].add(record.POS)
            except KeyError:
                data['reference'][str(record.REF)] = set([record.POS])

        logging.info("Completed %i out of %i vcfs", cnt, nof_vcfs)

    # check that there is no conflicting ref bases by verifying that the
    # intersection between two ref bases sets of positions is always empty
    for (contig, data) in all_contig_data.iteritems():
        ref_bases = data['reference'].keys()
        for i in range(0, len(ref_bases)):
            for j in range(0, len(ref_bases)):
                if i < j:
                    x = ref_bases[i]
                    y = ref_bases[j]
                    itscn = data['reference'][x] & data['reference'][y]
                    assert len(itscn) == 0, "FATAL ERROR: Two different ref bases for the same position: %s on contig %s" % (str(itscn), contig)

    """
    all_contig_data looks like this:
    {'gi|194097589|ref|NC_011035.1|': {'211696_H14256028001': {'-': set([12673,
                                                                         12674,
                                                                         12675,
                                                                         12676,
                                                                         12677]),
                                                                'A': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                'C': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                 ...
                                                               }
                                       }
                                       {'211697_H14333548001': {'-': set([12673,
                                                                         12674,
                                                                         12675,
                                                                         12676,
                                                                         12677]),
                                                                'A': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                'C': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                 ...
                                        ...

                                                               }
        ...
    }
    """

    # insert any filtering here ...
    logging.info("Filtering alignment.")

    if args["exclude"] or args["include"]:
        process_bed_file(args, all_contig_data)

    if args['remove_invariant_npos'] == True:
        for (contig, data) in all_contig_data.iteritems():
            # get all positions that are N and all others in all samples except the reference
            n_pos = set()
            var_pos = set()
            for sam in data.keys():
                if sam != 'reference':
                    for nuc in data[sam].keys():
                        if nuc == 'N':
                            n_pos.update(data[sam][nuc])
                        else:
                            var_pos.update(data[sam][nuc])
            # get positions that are onlys in nothing else in any sample
            only_n_pos = n_pos.difference(var_pos)
            # remove those positions from all samples including the reference
            for sam in data.keys():
                for nuc in data[sam].keys():
                    data[sam][nuc].difference_update(only_n_pos)

    if args['sample_Ns']:
        remove_samples(args, 'sample_Ns', 'N', all_contig_data)

    if args['sample_gaps']:
        remove_samples(args, 'sample_gaps', '-', all_contig_data)

    if args['column_Ns']:
        remove_columns(args['column_Ns'], 'N', all_contig_data)

    if args['column_gaps']:
        remove_columns(args['column_gaps'], '-', all_contig_data)

    # finished filtering
    logging.info("Filtering complete.")

    # output stats if required
    if args["with_stats"] is not None:
        logging.info("Calculating per position stats.")
        output_per_position_stats(args["with_stats"], all_contig_data, nof_vcfs)

    # output per sample stats
    logging.info("Writing per sample stats.")
    output_per_sample_stats(all_contig_data)

    # output now
    dSeqs = {}
    for (contig, data) in all_contig_data.iteritems():
        dAlign = {}
        # get all positions
        if args["reference"]:
            # all positions for whole contig when reference is required
            all_pos = {i + 1: i for i in range(len(args["reference"][contig]))}
        else:
            # get all positions for variant positions only else
            all_pos = set()
            for nuc in data['reference'].keys():
                all_pos.update(data['reference'][nuc])
            # all_pos is now a list ... and then a dict
            all_pos = sorted(all_pos)
            all_pos = {all_pos[i]: i for i in range(len(all_pos))}

        # 'initialise' sequence
        for sample_name in data.keys():
            if args["reference"]:
                # this is effectively a deepcopy, otherwise we're writing on the same
                # copy of the reference string for all samples
                # dAlign[sample_name] = args["reference"][contig]
                dAlign[sample_name] = []
                for s in args["reference"][contig]:
                    dAlign[sample_name].append(s)
            else:
                # initialies with 0's
                dAlign[sample_name] = ['0'] * len(all_pos)
                # set all bases to reference
                for nuc in data['reference'].keys():
                    for i in data['reference'][nuc]:
                        seq_pos = all_pos[i]
                        dAlign[sample_name][seq_pos] = nuc

            # overwrite reference positions where necessary
            for nuc in data[sample_name].keys():
                for i in data[sample_name][nuc]:
                    seq_pos = all_pos[i]
                    dAlign[sample_name][seq_pos] = nuc

            seq = ''.join(dAlign[sample_name])
            try:
                dSeqs[sample_name] += seq
            except KeyError:
                dSeqs[sample_name] = seq

    # write to file
    with open(args["out"], "w") as fp:
        # write seqs to file
        for name, seq in dSeqs.iteritems():
            fp.write(">%s\n%s\n" % (name, seq))

    return 0

# end of main --------------------------------------------------------------------------------------

def output_per_position_stats(filename, all_contig_data, nof_vcfs):
    """
    Output stats per position in the aligment.

    Parameters:
    -----------
    filename: str
        name of file to write to
    all_contig_data: dict
        all data as described abobe in main
    nof_vcfs: int
        number of samples excluding reference

    Returns:
    --------
    0
    """

    with open(filename, 'wb') as f:
        f.write("contig,position,mutations,n_frac,n_gaps\n")
        for contig, con_data in all_contig_data.iteritems():
            all_contig_pos = set()
            for sam in con_data.keys():
                for x in con_data[sam].values():
                    all_contig_pos.update(x)
            for pos in sorted(all_contig_pos):
                mut = 0.0
                ns = 0.0
                gaps = 0.0
                for sam in con_data.keys():
                    if sam == 'reference': continue
                    for m in ['A', 'C', 'G', 'T']:
                        try:
                            mut += 1.0 if pos in con_data[sam][m] else 0.0
                        except KeyError:
                            pass
                    try:
                        ns += 1.0 if pos in con_data[sam]['N'] else 0.0
                    except KeyError:
                        pass
                    try:
                        gaps += 1.0 if pos in con_data[sam]['-'] else 0.0
                    except KeyError:
                        pass
                f.write("%s,%i,%0.5f,%0.5f,%0.5f\n" % (contig,
                                                       pos,
                                                       mut/nof_vcfs,
                                                       ns/nof_vcfs,
                                                       gaps/nof_vcfs))

    return 0

# -------------------------------------------------------------------------------------------------

def output_per_sample_stats(all_contig_data):
    """
    Output stats per sample in the aligment.

    Parameters:
    -----------
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0
    """

    all_contigs = all_contig_data.keys()
    all_samples = set()
    for contig in all_contigs:
        all_samples.update(all_contig_data[contig].keys())
    for sam in all_samples:
        tot = 0
        ns = 0
        gaps = 0
        mut = 0
        mix = 0
        for contig in all_contigs:
            for k in all_contig_data[contig][sam].keys():
                tot += len(all_contig_data[contig][sam][k])
            try:
                ns += len(all_contig_data[contig][sam]['N'])
            except KeyError:
                pass
            try:
                gaps += len(all_contig_data[contig][sam]['-'])
            except KeyError:
                pass
            for k in ['A', 'C', 'G', 'T']:
                try:
                    mut += len(all_contig_data[contig][sam][k])
                except KeyError:
                    pass
            mixbases = set(all_contig_data[contig][sam].keys()).difference(['A', 'C', 'G', 'T', 'N', '-'])
            for m in mixbases:
                mix += len(all_contig_data[contig][sam][m])
        if sam == 'reference':
            mut = 0
        sys.stdout.write("%s\tN: %i, mut: %i, mix: %i, gap: %i, total: %i\n" % \
                         (sam, ns, mut, mix, gaps, tot))

    return 0

# --------------------------------------------------------------------------------------------------

def remove_columns(t, character, all_contig_data):
    """
    Remove columns form the alignment if the character [N, gap] is above
    the fraction t in a column.

    Parameters:
    -----------
    t: float
        threshold between 0.0 and 1.0
    character: str
        'N' or '-'
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0, but all_contig_data is updated
    """

    for (_, data) in all_contig_data.iteritems():
        all_pos = set()
        for nuc in data['reference']:
            all_pos.update(data['reference'][nuc])
        # number of samples not considering the reference
        nof_samples = len(data.keys()) - 1
        to_remove = set()
        for pos in all_pos:
            ns = 0
            for sam in data.keys():
                if sam != 'reference':
                    try:
                        if pos in data[sam][character]:
                            ns += 1
                    except KeyError:
                        # happens when sam doesn't have any Ns
                        pass
            if float(ns) / nof_samples > t:
                to_remove.add(pos)

        # remove postions
        if len(to_remove) > 0:
            for sam in data.keys():
                for nuc in data[sam].keys():
                    data[sam][nuc].difference_update(to_remove)

    return 0

# --------------------------------------------------------------------------------------------------

def process_bed_file(args, all_contig_data):
    """
    Process a bed file with position intervals to include or exclude.

    Parameters:
    -----------
    args: dict
        arguments dictionary
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0, but all_contig_data is updated
    """

    all_samples = set()
    for contig in all_contig_data.keys():
        all_samples.update(all_contig_data[contig].keys())

    incl = True if args['include'] else False

    with open(args["exclude"]) as fp:
        for line in fp:
            data = line.strip().split("\t")
            try:
                for sam in all_samples:
                    for nuc in all_contig_data[data[0]][sam].keys():
                        if incl == True:
                            all_contig_data[data[0]][sam][nuc].intersection_update(range(int(data[1]),
                                                                                         int(data[2]) + 1))
                        else:
                            all_contig_data[data[0]][sam][nuc].difference_update(range(int(data[1]),
                                                                                       int(data[2]) + 1))
            except KeyError:
                logging.error("Wrong contig in bed file: %s. Ignoring.", data[0])

    return 0

# --------------------------------------------------------------------------------------------------

def remove_samples(args, option, character, all_contig_data):
    """
    Remove samples form the alignment if the character [N, gap] is above
    a given fraction in in samples (relative to the size of the genome.)

    Parameters:
    -----------
    args: dict
        arguments dictionary
    option: str
        'sample_Ns' or 'sample_gaps'
    character: str
        'N' or '-'
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0, but all_contig_data is updated
    """

    # get the alignment length
    align_len = args["reflength"]

    # sum up the number of Ns or gaps for each sample
    ns_per_sample = {}
    for (_, data) in all_contig_data.iteritems():
        for sam in data.keys():
            if ns_per_sample.has_key(sam) == False:
                ns_per_sample[sam] = 0
            try:
                ns_per_sample[sam] += len(data[sam][character])
            except KeyError:
                # this happens when sam has no Ns
                pass

    # calculate proportion of Ns or gaps
    for sam in ns_per_sample.keys():
        ns_per_sample[sam] = float(ns_per_sample[sam]) / align_len

    # if this is set to auto calculate the threshold from the mean + 2 times the stddev
    t = 0.0
    if args[option] == 'auto':
        m = sum(ns_per_sample.values()) / len(ns_per_sample)
        ssd = sum((x-m)**2 for x in ns_per_sample.values())
        variance = ssd / len(ns_per_sample)
        sd = sqrt(variance)
        t = m + (args['sample_Ns_gaps_auto_factor']*sd)
        logging.info("Mean number of %ss per sample is %.6f and standard deviation is %.6f.",
                     'gap' if character == '-' else character,
                     m, sd)
    else:
        t = args[option]
    logging.info("Threshold for %ss per sample set to %.6f. Removing samples that have more.",
                 'gap' if character == '-' else character, t)

    # remove all samples that are bigger than the threshold
    removals = False
    for sam in ns_per_sample.keys():
        if ns_per_sample[sam] > t:
            for (_, data) in all_contig_data.iteritems():
                logging.info("Removing sample %s, because it has %.6f %ss",
                             sam,
                             ns_per_sample[sam],
                             'gap' if character == '-' else character)
                del data[sam]
                removals = True

    # tidy up
    if removals == True:
        for (contig, data) in all_contig_data.iteritems():
            # get all positions in the reference and all positions in all other samples
            ref_pos = set()
            var_pos = set()
            for sam in data.keys():
                if sam == 'reference':
                    for nuc in data[sam].keys():
                        ref_pos.update(data[sam][nuc])
                else:
                    for nuc in data[sam].keys():
                        var_pos.update(data[sam][nuc])
            # get all positions that are only in the reference
            only_ref_pos = ref_pos.difference(var_pos)
            # remove the positions only in the reference from the reference
            # because we probably removed the sample with a variant at those positions
            for nuc in data['reference'].keys():
                data['reference'][nuc].difference_update(only_ref_pos)

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    exit(main(vars(get_args().parse_args())))
