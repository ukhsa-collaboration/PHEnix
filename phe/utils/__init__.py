'''
:Date: 12May2016
:Author: Publich Health England
'''

import vcf
from bintrees import FastRBTree
import math
import logging

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
        return "N: %i, mut: %i, mix: %i, gap: %i, total: %i" % \
               (self.N, self.mut, self.mix, self.gap, self.total)

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

    # First pass to get the references and the positions to be analysed.
    for vcf_in in dArgs['input']:

        reader = vcf.Reader(filename=vcf_in)

        # Get the sample name from the VCF file (usually the read group).
        sample_name = reader.samples[0]
        aSampleNames.append(sample_name)

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
                # mix = get_mixture(record, args.with_mixtures)
                # Currently we are only using first filter to call consensus.
                extended_code = "N"
                if extended_code == "N":
                    position_data["stats"].N += 1
                position_data[sample_name] = extended_code

            else:
                # filter fail; code as N for consistency
                position_data[sample_name] = "N"
                position_data["stats"].N += 1

    # finished parsing

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

    # initialise empty matrix
    dist_mat = {}
    for i, sample_1 in enumerate(aSampleNames):
        dist_mat[sample_1] = {}
        for j, sample_2 in enumerate(aSampleNames):
            if j < i:
                continue
            if dArgs['substitution'] == 'k80' or dArgs['substitution'] == 't93':
                dist_mat[sample_1][sample_2] = [0.0, 0.0]
            elif dArgs['substitution'] == 'tn84':
                dist_mat[sample_1][sample_2] = {'A': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                                                'C': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                                                'G': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                                                'T': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}}
            else:
                dist_mat[sample_1][sample_2] = 0.0

    for sContig in avail_pos.keys():
        for pos in avail_pos[sContig]:
            ref_base = avail_pos[sContig][pos].get("reference")
            for i, sample_1 in enumerate(aSampleNames):
                s1_base = avail_pos[sContig][pos].get(sample_1, ref_base)
                if dValChars.get(s1_base.upper(), None) == None:
                    continue
                for j, sample_2 in enumerate(aSampleNames):
                    if j <= i:
                        continue
                    s2_base = avail_pos[sContig][pos].get(sample_2, ref_base)
                    if dValChars.get(s2_base.upper(), None) == None:
                        continue


                    if dArgs['substitution'] == 'k80' or dArgs['substitution'] == 't93':
                        k = get_difference_value(s1_base.upper(), s2_base.upper(), dArgs['substitution'])
                        dist_mat[sample_1][sample_2][0] += k[0]
                        dist_mat[sample_1][sample_2][1] += k[1]
                    elif dArgs['substitution'] == 'tn84':
                        if s1_base > s2_base:
                            s1_base, s2_base = s2_base, s1_base
                        dist_mat[sample_1][sample_2][s1_base][s2_base] += 1.0
                    else:
                        k = get_difference_value(s1_base.upper(), s2_base.upper(), dArgs['substitution'])
                        dist_mat[sample_1][sample_2] += k


    if dArgs['substitution'] == 'jc69':
        dist_mat = normalise_jc69(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 'k80':
        dist_mat = normalise_k80(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 'tn84':
        dist_mat = normalise_tn84(dist_mat, dArgs['refgenome'], aSampleNames)
    elif dArgs['substitution'] == 't93':
        dist_mat = normalise_t93(dist_mat, dArgs['refgenome'], aSampleNames)
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
                Q = d[sample_2][sample_1][1] / flGenLen
                p = sum(d[sample_2][sample_1]) / flGenLen
                x1 = -h * math.log(1.0 - p/h - Q)
                x2 = 0.5 * (1.0 - h) * math.log(1.0 - (2.0 * Q))
                d[sample_2][sample_1] = x1 - x2
            elif i == j:
                d[sample_2][sample_1] = 0.0

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
                iNofDiff = getTotalNofDiff_tn84(d[sample_2][sample_1])
                p = iNofDiff / flGenLen
                sum2 = 0.0
                for m, nuc1 in enumerate(['A', 'C', 'G', 'T']):
                    for n, nuc2 in enumerate(['A', 'C', 'G', 'T']):
                        if m < n:
                            flFreqNucPair = d[sample_2][sample_1][nuc1][nuc2] / flGenLen
                            sum2 += (flFreqNucPair ** 2.0) / (2.0 * dRefFreq[nuc1] * dRefFreq[nuc2])
                b = 0.5 * ((1.0 - sum1) + (p * p / sum2))
                dist = -b * math.log(1.0 - (p / b))
                d[sample_2][sample_1] = dist
            elif i == j:
                d[sample_2][sample_1] = 0.0
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
                P = d[sample_2][sample_1][0] / flGenLen
                Q = d[sample_2][sample_1][1] / flGenLen
                w1 = 1.0 - (2.0 * P) - Q
                w2 = 1.0 - (2.0 * Q)
                x = (-0.5 * math.log(w1)) - (0.25 * math.log(w2))
                d[sample_2][sample_1] = x
            elif i == j:
                d[sample_2][sample_1] = 0.0

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
                p = d[sample_2][sample_1] / flGenLen
                x = (-3.0/4.0) * math.log(1.0-((4.0/3.0) * p))
                d[sample_2][sample_1] = x

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
