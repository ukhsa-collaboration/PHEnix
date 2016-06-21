'''
:Date: 12May2016
:Author: Public Health England
'''

import sys
import os
import vcf
from bintrees import FastRBTree
import math
import logging
from numpy import std
from numpy import mean

from pprint import pprint
from collections import OrderedDict

dValChars = {'A': 1, 'C': 1, 'G': 1, 'T': 1}

# --------------------------------------------------------------------------------------------------

# G -> A = transition = (1, 0)
# A -> T = transversion = (0, 1)
dK80 = {'A': {'A': [0.0, 0.0], 'C': [0.0, 1.0], 'G': [1.0, 0.0], 'T': [0.0, 1.0]},
        'C': {'A': [0.0, 1.0], 'C': [0.0, 0.0], 'G': [0.0, 1.0], 'T': [1.0, 0.0]},
        'G': {'A': [1.0, 0.0], 'C': [0.0, 1.0], 'G': [0.0, 0.0], 'T': [0.0, 1.0]},
        'T': {'A': [0.0, 1.0], 'C': [1.0, 0.0], 'G': [0.0, 1.0], 'T': [0.0, 0.0]}}

# --------------------------------------------------------------------------------------------------

class base_stats(object):
    def __init__(self):
        self.N = 0
        self.mut = 0
        self.gap = 0
        self.mix = 0
        self.total = 0

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


# --------------------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------------------

def precompute_snp_densities(avail_pos, sample_names, ref, k):
    """
    avail_pos:
    {'gi|194097589|ref|NC_011035.1|':
        FastRBTree({2329: {'stats': <vcf2distancematrix.base_stats object at 0x40fb590>,
                           'reference': 'A',
                           '211700_H15498026501': 'C',
                           '211701_H15510030401': 'C',
                           '211702_H15522021601': 'C'},
                    3837: {'211700_H15498026501': 'G',
                           'stats': <vcf2distancematrix.base_stats object at 0x40fbf90>,
                           '211701_H15510030401': 'G',
                           'reference': 'T',
                           '211702_H15522021601': 'G'},
                    4140: {'211700_H15498026501': 'A',
                           'stats': <vcf2distancematrix.base_stats object at 0x40fb790>,
                           '211701_H15510030401': 'A',
                           'reference': 'G',
                           '211702_H15522021601': 'A'}})}

    returns dDen:
    {'gi|194097589|ref|NC_011035.1|': {2329: {'211700_H15498026501': 0.001,
                                              '211701_H15510030401': 0.001,
                                              '211702_H15522021601': 0.001},
                                       3837: {'211700_H15498026501': 0.002,
                                              '211701_H15510030401': 0.002,
                                              '211702_H15522021601': 0.002},
                                       4140: {'211700_H15498026501': 0.002,
                                              '211701_H15510030401': 0.002,
                                              '211702_H15522021601': 0.002},
                                       12159: {'211700_H15498026501': 0.009,
                                               '211701_H15510030401': 0.008,
                                               '211702_H15522021601': 0.008},
                                       12309: {'211700_H15498026501': 0.149,
                                               '211701_H15510030401': 0.148,
                                               '211702_H15522021601': 0.148}},
     'per_samplepair_threshold': {'211700_H15498026501': {},
                                  '211701_H15510030401': {'211700_H15498026501': 0.20127650534341832},
                                  '211702_H15522021601': {'211700_H15498026501': 0.20836662410991291,
                                                          '211701_H15510030401': 0.20455899117993009}}}
    """

    iWINSIZE = 1000
    flWINSIZE = float(iWINSIZE)

    dSNPWins = {}
    (dRefFreq, flGenLen) = get_ref_freqs(ref)
    for i, sname1 in enumerate(sample_names):
        dSNPWins[sname1] = []
        for j in range(0, int(flGenLen), iWINSIZE):
            dSNPWins[sname1].append(0.0)

    for contig, oBT in avail_pos.items():
        for iPos in oBT:
            ref_base = oBT[iPos]['reference']
            if dValChars.get(ref_base, None) == None:
                continue
            for sname in sample_names:
                try:
                    sam_base = oBT[iPos][sname]
                except KeyError:
                    continue
                if dValChars.get(sam_base, None) == None:
                    continue
                if sam_base != ref_base:
                    dSNPWins[sname][iPos/iWINSIZE] += 1/flWINSIZE

    dDen = {}
    dPerSample = {}
    for contig, oBT in avail_pos.items():
        dDen[contig] = {}
        for iPos in oBT:
            ref_base = oBT[iPos]['reference']
            if dValChars.get(ref_base, None) == None:
                continue
            for sname in sample_names:
                try:
                    sam_base = oBT[iPos][sname]
                except KeyError:
                    continue
                if dValChars.get(sam_base, None) == None:
                    continue
                if sam_base != ref_base:
                    flSNPsInWin = 0.0
                    for x in range((iPos - ((iWINSIZE/2)-1)), (iPos + (iWINSIZE/2) + 1)):
                        try:
                            window_base = oBT[x][sname]
                            if dValChars.get(window_base, None) == None:
                                continue
                            window_ref = oBT[x]['reference']
                            if dValChars.get(window_ref, None) == None:
                                continue
                            if window_base != window_ref:
                                flSNPsInWin += 1.0
                        except KeyError:
                            continue
                    flDnsty = flSNPsInWin / flWINSIZE
                    try:
                        dDen[contig][iPos][sname] = flDnsty
                    except KeyError:
                        dDen[contig][iPos] = {}
                        dDen[contig][iPos][sname] = flDnsty

    dDen['per_samplepair_threshold'] = {}
    for i, sname1 in enumerate(sample_names):
        dDen['per_samplepair_threshold'][sname1] = {}
        for j, sname2 in enumerate(sample_names):
            if j < i:
                density_windows = []
                if sname1 != 'reference':
                    density_windows += dSNPWins[sname1]
                if sname2 != 'reference':
                    density_windows += dSNPWins[sname2]
                tmp_thresh = mean(density_windows) + (k * std(density_windows))
                dDen['per_samplepair_threshold'][sname1][sname2] = max(1.0/flWINSIZE, tmp_thresh)

    return dDen

# --------------------------------------------------------------------------------------------------

def parse_vcf_files(dArgs, avail_pos, aSampleNames):
    '''
    Function
    Parameters
    ----------

    Returns
    -------

    '''

    empty_tree = FastRBTree()

    exclude = {}
    include = {}

    if dArgs["exclude"] != None or dArgs["include"] != None:
        pos = {}
        chr_pos = []
        bed_file = dArgs["include"] if dArgs["include"] != None else dArgs["exclude"]

        with open(bed_file) as fp:
            for line in fp:
                data = line.strip().split("\t")
                chr_pos += [(i, False,) for i in xrange(int(data[1]), int(data[2]) + 1)]
                try:
                    pos[data[0]] += chr_pos
                except KeyError:
                    pos[data[0]] = chr_pos
        pos = {chrom: FastRBTree(l) for chrom, l in pos.items()}
        if dArgs["include"]:
            include = pos
        else:
            exclude = pos

    dSamNms = {'reference': None}

    # First pass to get the references and the positions to be analysed.
    for vcf_in in dArgs['input']:

        reader = vcf.Reader(filename=vcf_in)

        # Get the sample name from the VCF file (usually the read group).
        sample_name = reader.samples[0]
        assert dSamNms.has_key(sample_name) == False, "ERROR: %s is not a unique sample name in the set of vcfs" % (sample_name)
        dSamNms[sample_name] = None

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
            # if include.get(record.CHROM, empty_tree).get(record.POS, False) or \
            #    not exclude.get(record.CHROM, empty_tree).get(record.POS, True):
                continue

            try:
                _ = avail_pos[record.CHROM]
            except KeyError:
                avail_pos[record.CHROM] = FastRBTree()

            # Setup the position data to contain reference and stats.
            if avail_pos[record.CHROM].get(record.POS, None) is None:
                avail_pos[record.CHROM].insert(record.POS, {"reference": str(record.REF),
                                                            "stats": base_stats()})

            position_data = avail_pos[record.CHROM].get(record.POS)

            assert len(position_data["reference"]) == 1, "Reference base must be singluar: in %s found %s @ %s" % (position_data["reference"], sample_name, record.POS)

            # IF this is uncallable genotype, add gap "-"
            if record.is_uncallable:
                # TODO: Mentioned in issue: #7(gitlab)
                position_data[sample_name] = "-"
                # Update stats
                position_data["stats"].gap += 1

            elif not record.FILTER:
                # If filter PASSED!
                # Make sure the reference base is the same. Maybe a vcf from different species snuck in here?!
                assert str(record.REF) == position_data["reference"] or str(record.REF) == 'N' or position_data["reference"] == 'N', "SOMETHING IS REALLY WRONG because reference for the same position is DIFFERENT! %s in %s (%s, %s)" % (record.POS, vcf_in, str(record.REF), position_data["reference"])
                # update position_data['reference'] to a real base if possible
                if position_data['reference'] == 'N' and str(record.REF) != 'N':
                    position_data['reference'] = str(record.REF).upper()
                if record.is_snp:
                    if len(record.ALT) > 1:
                        logging.info("POS %s passed filters but has multiple alleles REF: %s, ALT: %s. Inserting N" % (record.POS, str(record.REF), str(record.ALT)))
                        position_data[sample_name] = "N"
                        position_data["stats"].N += 1

                    else:
                        position_data[sample_name] = str(record.ALT[0]).upper()
                        position_data["stats"].mut += 1

            # Filter(s) failed
            elif record.is_snp:
                position_data[sample_name] = 'N'
                position_data["stats"].N += 1
            else:
                # filter fail; code as N for consistency
                position_data[sample_name] = "N"
                position_data["stats"].N += 1
    # finished parsing

    for sn in dSamNms.keys():
        aSampleNames.append(sn)

    return 0

# --------------------------------------------------------------------------------------------------

def get_dist_mat(aSampleNames, avail_pos, dArgs):
    """
    avail_pos:
    {'gi|194097589|ref|NC_011035.1|':
        FastRBTree({2329: {'stats': <vcf2distancematrix.base_stats object at 0x40fb590>,
                           'reference': 'A',
                           '211700_H15498026501': 'C',
                           '211701_H15510030401': 'C',
                           '211702_H15522021601': 'C'},
                    3837: {'211700_H15498026501': 'G',
                           'stats': <vcf2distancematrix.base_stats object at 0x40fbf90>,
                           '211701_H15510030401': 'G',
                           'reference': 'T',
                           '211702_H15522021601': 'G'},
                    4140: {'211700_H15498026501': 'A',
                           'stats': <vcf2distancematrix.base_stats object at 0x40fb790>,
                           '211701_H15510030401': 'A',
                           'reference': 'G',
                           '211702_H15522021601': 'A'}})}
    """

    dDen = None
    dRemovals = None
    if dArgs['remove_recombination'] == True:
        dDen = precompute_snp_densities(avail_pos, aSampleNames, dArgs['refgenome'], dArgs['k'])
        dRemovals = {}
        for i, sample_1 in enumerate(aSampleNames):
            dRemovals[sample_1] = {}
            for j, sample_2 in enumerate(aSampleNames):
                if j <= i:
                    dRemovals[sample_1][sample_2] = 0


    # initialise empty matrix
    dist_mat = {}
    for i, sample_1 in enumerate(aSampleNames):
        dist_mat[sample_1] = {}
        for j, sample_2 in enumerate(aSampleNames):
            if j <= i:
                if dArgs['substitution'] == 'k80' or dArgs['substitution'] == 't93':
                    dist_mat[sample_1][sample_2] = [0.0, 0.0]
                elif dArgs['substitution'] == 'tn84':
                    dist_mat[sample_1][sample_2] = {'A': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                                                    'C': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                                                    'G': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                                                    'T': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}}
                else:
                    dist_mat[sample_1][sample_2] = 0.0
            else: # j > i
                pass

    aStats = []
    for sContig in avail_pos.keys():
        for pos in avail_pos[sContig]:
            ref_base = avail_pos[sContig][pos].get("reference")
            for i, sample_1 in enumerate(aSampleNames):
                s1_base = avail_pos[sContig][pos].get(sample_1, ref_base)
                # consider only differences between valid characters -> pairwise deletion
                if dValChars.get(s1_base.upper(), None) == None:
                    continue
                for j, sample_2 in enumerate(aSampleNames):
                    if j < i:
                        s2_base = avail_pos[sContig][pos].get(sample_2, ref_base)
                        # consider only differences between valid characters -> pairwise deletion
                        if dValChars.get(s2_base.upper(), None) == None:
                            continue

                        if dDen != None and s1_base != s2_base:

                            flLocalDensity1 = dDen[sContig][pos].get(sample_1, 0.0)
                            flLocalDensity2 = dDen[sContig][pos].get(sample_2, 0.0)
                            flPairThreshold = dDen['per_samplepair_threshold'][sample_1][sample_2]
                            # aStats.append("%s\t%s\t%i\t%.3f\t%.3f\t%.3f\n" % (sample_1, sample_2, pos, flLocalDensity1, flLocalDensity2, flPairThreshold))
                            if (flLocalDensity1 > flPairThreshold) or (flLocalDensity2 > flPairThreshold):
                                aStats.append("%s\t%s\t%i\t%.3f\t%.3f\t%.3f\n" % (sample_1, sample_2, pos, flLocalDensity1, flLocalDensity2, flPairThreshold))
                                dRemovals[sample_1][sample_2] += 1
                                continue

                        if dArgs['substitution'] == 'k80' or dArgs['substitution'] == 't93':
                            k = get_difference_value(s1_base.upper(), s2_base.upper(), dArgs['substitution'])
                            dist_mat[sample_1][sample_2][0] += k[0]
                            dist_mat[sample_1][sample_2][1] += k[1]
                        elif dArgs['substitution'] == 'tn84':

                            if s1_base > s2_base:
                                # don't do this:
                                # s1_base, s2_base = s2_base, s1_base
                                # dist_mat[sample_1][sample_2][s1_base][s2_base] += 1.0
                                # => messes up the references
                                dist_mat[sample_1][sample_2][s2_base][s1_base] += 1.0
                            else:
                                dist_mat[sample_1][sample_2][s1_base][s2_base] += 1.0

                        elif dArgs['substitution'] == 'number_of_differences' or dArgs['substitution'] == 'jc69':

                            k = get_difference_value(s1_base.upper(), s2_base.upper(), dArgs['substitution'])
                            dist_mat[sample_1][sample_2] += k
                        else:
                            raise NotImplementedError
                    else: # j >= i
                        pass

    # write additional stats if required
    if dArgs['remove_recombination'] == True and dArgs['with_stats'] == True:
        sOutBase = os.path.splitext(dArgs['out'])[0]
        with open("%s.removals.tsv" % (sOutBase), 'w') as fOut:
            for i, sample_1 in enumerate(aSampleNames):
                row = sample_1
                for j, sample_2 in enumerate(aSampleNames):
                    if j < i:
                        row += "%s%i" % ('\t', dRemovals[sample_1][sample_2])
                fOut.write("%s\n" % row)
        with open("%s.proportion_removed.tsv" % (sOutBase), 'w') as fOut:
            for i, sample_1 in enumerate(aSampleNames):
                row = sample_1
                for j, sample_2 in enumerate(aSampleNames):
                    if j < i:
                        try:
                            row += "%s%f" % ('\t', dRemovals[sample_1][sample_2] \
                                                 / (dist_mat[sample_1][sample_2] \
                                                 + dRemovals[sample_1][sample_2]))
                        except ZeroDivisionError:
                            row += "\tNAN"
                fOut.write("%s\n" % row)
        with open("%s.removed_snps.tsv" % (sOutBase), 'w') as fOut:
            fOut.write("sample_1\tsample_2\tposition\tdensity_in_sample1\tdensity_in_sample2\tthreshold_for_pair\n")
            for sLine in aStats:
                fOut.write(sLine)

    if dArgs['substitution'] == 'jc69':
        dist_mat = normalise_jc69(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 'k80':
        dist_mat = normalise_k80(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 'tn84':
        dist_mat = normalise_tn84(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 't93':
        dist_mat = normalise_t93(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 'number_of_differences':
        pass
    else:
        raise NotImplementedError

    return dist_mat

# --------------------------------------------------------------------------------------------------

def normalise_t93(d, ref, names):
    """
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 16
    """

    (dRefFreq, flGenLen) = get_ref_freqs(ref)
    flGC = dRefFreq['C'] + dRefFreq['G']
    h = 2.0 * flGC * (1.0 - flGC)

    for i, sample_1 in enumerate(names):
        for j, sample_2 in enumerate(names):
            if j < i:
                Q = d[sample_1][sample_2][1] / flGenLen
                p = sum(d[sample_1][sample_2]) / flGenLen
                x1 = -h * math.log(1.0 - p/h - Q)
                x2 = 0.5 * (1.0 - h) * math.log(1.0 - (2.0 * Q))
                d[sample_1][sample_2] = x1 - x2
            elif i == j:
                d[sample_1][sample_2] = 0.0

    return d

# --------------------------------------------------------------------------------------------------

def normalise_tn84(d, ref, names):
    """
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 13 and 14

    d = {'211701_H15510030401': {'211700_H15498026501': {'A': {'A': 1152.0,
                                                               'C': 114.0,
                                                               'G': 545.0,
                                                               'T': 35.0},
                                                         'C': {'A': 0.0,
                                                               'C': 1233.0,
                                                               'G': 108.0,
                                                               'T': 467.0},
                                                         'G': {'A': 0.0,
                                                               'C': 0.0,
                                                               'G': 1283.0,
                                                               'T': 100.0},
                                                         'T': {'A': 0.0,
                                                               'C': 0.0,
                                                               'G': 0.0,
                                                               'T': 1177.0}}, ...}, ...}
    """

    (dRefFreq, flGenLen) = get_ref_freqs(ref)

    sum1 = sum([x*x for x in dRefFreq.values()])

    for i, sample_1 in enumerate(names):
        for j, sample_2 in enumerate(names):
            if j < i:
                iNofDiff = getTotalNofDiff_tn84(d[sample_1][sample_2])

                p = iNofDiff / flGenLen
                sum2 = 0.0
                for m, nuc1 in enumerate(['A', 'C', 'G', 'T']):
                    for n, nuc2 in enumerate(['A', 'C', 'G', 'T']):
                        if m < n:
                            flFreqNucPair = d[sample_1][sample_2][nuc1][nuc2] / flGenLen
                            sum2 += ((flFreqNucPair ** 2.0) / (2.0 * dRefFreq[nuc1] * dRefFreq[nuc2]))
                if sum2 == 0.0: # samples are "identical"
                    d[sample_1][sample_2] = 0.0
                else:
                    b = 0.5 * ((1.0 - sum1) + (p * p / sum2))
                    d[sample_1][sample_2] = -b * math.log(1.0 - (p / b))
            elif i == j:
                d[sample_1][sample_2] = 0.0
            else: # j > i
                pass

    return d

# --------------------------------------------------------------------------------------------------

def normalise_k80(d, ref, names):
    """
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 9
    """

    (_, flGenLen) = get_ref_freqs(ref)

    for i, sample_1 in enumerate(names):
        for j, sample_2 in enumerate(names):
            if j < i:
                P = d[sample_1][sample_2][0] / flGenLen
                Q = d[sample_1][sample_2][1] / flGenLen
                w1 = 1.0 - (2.0 * P) - Q
                w2 = 1.0 - (2.0 * Q)
                d[sample_1][sample_2] = (-0.5 * math.log(w1)) - (0.25 * math.log(w2))
            elif i == j:
                d[sample_1][sample_2] = 0.0

    return d

# --------------------------------------------------------------------------------------------------

def normalise_jc69(d, ref, names):
    """
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 7
    """

    flGenLen = 0.0
    with open(ref, 'r') as fRef:
        for sLine in fRef:
            if sLine.startswith(">"):
                continue
            else:
                flGenLen += len(sLine.strip())

    logging.info("genome length: %s" % flGenLen)

    for i, sample_1 in enumerate(names):
        for j, sample_2 in enumerate(names):
            if j < i:
                p = d[sample_1][sample_2] / flGenLen
                d[sample_1][sample_2] = (-3.0/4.0) * math.log(1.0-((4.0/3.0) * p))

    return d

# --------------------------------------------------------------------------------------------------

def get_difference_value(s1_base, s2_base, sSubs):
    """
    todo: docstring
    """
    if sSubs == 'number_of_differences' or sSubs == 'jc69':

        try:
            _ = dValChars[s1_base]
        except KeyError:
            return 0.0
        try:
            _ = dValChars[s2_base]
        except KeyError:
            return 0.0

        if s1_base == s2_base:
            return 0.0
        else:
            return 1.0
    elif sSubs == 'k80' or sSubs == 't93':

        try:
            _ = dValChars[s1_base]
        except KeyError:
            return [0.0, 0.0]
        try:
            _ = dValChars[s2_base]
        except KeyError:
            return [0.0, 0.0]

        return dK80[s1_base][s2_base]
    else:
        raise NotImplementedError

# --------------------------------------------------------------------------------------------------

def getTotalNofDiff_tn84(d):
    """
    sum up total number of differences for a dict like this:
    {'A': {'A': 1152.0, 'C': 114.0, 'G': 545.0, 'T': 35.0},
     'C': {'A': 0.0, 'C': 1233.0, 'G': 108.0, 'T': 467.0},
     'G': {'A': 0.0, 'C': 0.0, 'G': 1283.0, 'T': 100.0},
     'T': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 1177.0}}
     = 1369.0 in this case
    """
    t = 0.0
    for k, v in d.items():
        t += sum([y if x != k else 0.0 for x, y in v.items()])
    return t

# --------------------------------------------------------------------------------------------------

def get_ref_freqs(ref):

    # get frequency of a, c, g, and t in ref
    dRefFreq = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
    flGenLen = 0.0
    with open(ref, 'r') as fRef:
        for sLine in fRef:
            if sLine.startswith(">"):
                continue
            else:
                sLine = sLine.upper().strip()
                flGenLen += len(sLine)
                for n in sLine:
                    try:
                        dRefFreq[n] += 1.0
                    except KeyError:
                        pass
    for (i, j) in dRefFreq.items():
        dRefFreq[i] = j / flGenLen

    logging.info("genome length: %s" % flGenLen)
    return (dRefFreq, flGenLen)

# --------------------------------------------------------------------------------------------------
