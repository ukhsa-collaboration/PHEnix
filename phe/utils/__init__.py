from collections import OrderedDict
import logging
from math import floor
import math
import multiprocessing
import os
from pprint import pprint
import sys

from bintrees import FastRBTree
from numpy import mean
from numpy import std
import vcf

HAVE_SCIPY = True
try:
    from scipy.stats import binom_test
except ImportError:
    HAVE_SCIPY = False

HAVE_PSUTIL = True
try:
    from psutil import virtual_memory
except ImportError:
    HAVE_PSUTIL = False

dValChars = {'A': 1, 'C': 1, 'G': 1, 'T': 1}

# --------------------------------------------------------------------------------------------------

# G -> A = transition = (1, 0)
# A -> T = transversion = (0, 1)
dK80 = {'A': {'A': [0.0, 0.0], 'C': [0.0, 1.0], 'G': [1.0, 0.0], 'T': [0.0, 1.0]},
        'C': {'A': [0.0, 1.0], 'C': [0.0, 0.0], 'G': [0.0, 1.0], 'T': [1.0, 0.0]},
        'G': {'A': [1.0, 0.0], 'C': [0.0, 1.0], 'G': [0.0, 0.0], 'T': [0.0, 1.0]},
        'T': {'A': [0.0, 1.0], 'C': [1.0, 0.0], 'G': [0.0, 1.0], 'T': [0.0, 0.0]}}

# --------------------------------------------------------------------------------------------------

class BaseStats(object):
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
    """Is the Record uncallable? Currently the record is **uncallable** iff:

     * GT field is **./.**
     * **LowQual** is in the filter.

    Returns
    -------
    uncall: bool
        True if any of the above items are true, False otherwise.
    """

    uncall = False
    try:
        if record.samples[0].data.GT in ("./.", None):
            uncall = True
    except AttributeError:
        uncall = None

    if record.FILTER is not None and "LowQual" in record.FILTER:
        uncall = True

    return uncall

# --------------------------------------------------------------------------------------------------

def precompute_snp_densities(avail_pos, sample_names, args):
    """
    Precompute the number of differences around each difference between each pair of samples

    Parameters
    ----------
    avail_pos: dict
        data structure that contains the information on all available positions, like this:
        {'gi|194097589|ref|NC_011035.1|':
        FastRBTree({2329: {'stats': <vcf2distancematrix.BaseStats object at 0x40fb590>,
                           'reference': 'A',
                           '211700_H15498026501': 'C',
                           '211701_H15510030401': 'C',
                           '211702_H15522021601': 'C'},
                    3837: {'211700_H15498026501': 'G',
                           'stats': <vcf2distancematrix.BaseStats object at 0x40fbf90>,
                           '211701_H15510030401': 'G',
                           'reference': 'T',
                           '211702_H15522021601': 'G'},
                    4140: {'211700_H15498026501': 'A',
                           'stats': <vcf2distancematrix.BaseStats object at 0x40fb790>,
                           '211701_H15510030401': 'A',
                           'reference': 'G',
                           '211702_H15522021601': 'A'}})}
    sample_names: list
        list of sample names
    args: dict
        input parameter dictionary as created by get_args()

    Returns
    -------
    dDen: dict
        contains the differences between a pair in a window of given size
        around each difference of the pair
        {'diffs': {'187534_H153520399-1': {'187534_H153520399-1': 0,
                                           '187536_H154060132-1': 1609,
                                           '189918_H154320283-2': 295,
                                           '205683_H15352039901': 0,
                                           '211698_H15464036401': 298,
                                           '211700_H15498026501': 298,
                                           '211701_H15510030401': 1621,
                                           '211702_H15522021601': 297,
                                           '211703_H15534021301': 1632,
                                           'reference': 4045},
                   '187536_H154060132-1': {'187536_H154060132-1': 0,
                                           '205683_H15352039901': 1605,
                                           '211698_H15464036401': 1353,
                                           '211701_H15510030401': 1,
                                           '211702_H15522021601': 1351,
                                           '211703_H15534021301': 7,
                                           'reference': 5041}
                                           ...,
                                           },
         'gi|194097589|ref|NC_011035.1|': {'187534_H153520399-1': {'187536_H154060132-1': {55959: 1,
                                                                                           56617: 1,
                                                                                           157165: 1,
                                                                                           279950: 3,
                                                                                           279957: 3,
                                                                                           279959: 3,
                                                                                           608494: 22,
                                                                                           608537: 23,
                                                                                           608551: 23,
                                                                                           608604: 23,
                                                                                           608617: 24,
                                                                                           ...,}
                                                                    '189918_H154320283-2': {27696: 1,
                                                                                            55959: 1,
                                                                                            56617: 2,
                                                                                            56695: 2,
                                                                                            279950: 3,
                                                                                            279957: 3,
                                                                                            279959: 3,
                                                                                            520610: 1,
                                                                                            608494: 22,
                                                                                            ...,
                                                                                           }}}}
    """

    (_, flGenLen) = get_ref_freqs(args['refgenome'], len_only=True)

    dDen = {}
    dDen['diffs'] = {}
    for i, sample_1 in enumerate(sample_names):
        dDen['diffs'][sample_1] = {}
        for j, sample_2 in enumerate(sample_names):
            if j <= i:
                dDen['diffs'][sample_1][sample_2] = 0

    pool = multiprocessing.Pool(args['threads'])

    parameters = []

    for contig, oBT in avail_pos.items():
        dDen[contig] = {}
        for i, sample_1 in enumerate(sample_names):
            dDen[contig][sample_1] = {}
            for j, sample_2 in enumerate(sample_names):
                if j < i:
                    parameters.append((sample_1, sample_2, oBT, flGenLen,))

    results = pool.map(_get_sample_pair_densities, parameters)

    # Close the pool and wait for all tasks to be completed.
    pool.close()
    pool.join()

    for contig, oBT in avail_pos.items():
        for i, sample_1 in enumerate(sample_names):
            for j, sample_2 in enumerate(sample_names):
                if j < i:
                    (x, y) = results.pop(0)
                    dDen['diffs'][sample_1][sample_2] += x
                    dDen[contig][sample_1][sample_2] = y

    # debug
    # sOutBase = os.path.splitext(args['out'])[0]
    # with open('%s_dDen.txt' % (sOutBase), 'w') as f:
    #     pprint(dDen, f)

    return dDen

# --------------------------------------------------------------------------------------------------

def _get_sample_pair_densities(ARGS):
    """
    This wrapper is required to call a funtion with pool.map that accepts >1 parameter.
    (see http://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
    answers by J.F. Sebastian and imotai)

    Parameters
    ----------
    ARGS: tuple
        tuple with parameters
    Returns
    -------
    call to get_sample_pair_densities with parameters unpacked
    """

    return get_sample_pair_densities(*ARGS)

# --------------------------------------------------------------------------------------------------

def get_sample_pair_densities(sample_1, sample_2, oBT, flGenLen):
    '''
    Function to calculate the differecnes in a window of a given size around is
    difference for a given pair

    Parameters
    ----------
    sample_1: str
        name of sample 1
    sample_2: str
        name of sample 2
    oBT: obj
        bintree object that contains all information for all available
        positions for a given contig
    flGenLen: float
        reference genome length
    Returns
    -------
    (diffs, d): tuple
        diffs: int
            total number of differences between the pair
        d: dict
            dict with position of difference as key and number of differences
            in window around it as value
    '''

    # first find out the total number of diff between the two samples
    diffs = 0
    for pos in oBT:
        ref_base = oBT[pos].get("reference")
        s1_base = oBT[pos].get(sample_1, ref_base)
        # consider only differences between valid characters -> pairwise deletion
        if dValChars.get(s1_base.upper(), None) == None:
            continue
        s2_base = oBT[pos].get(sample_2, ref_base)
        # consider only differences between valid characters -> pairwise deletion
        if dValChars.get(s2_base.upper(), None) == None:
            continue
        if s1_base != s2_base:
            diffs += 1


    d = {}

    # use win size of avg nof nucleotides per difference, i.e. each window should have one diff on average
    # for very similar sample pairs (very few diffs) do not impose max win size
    # if >1% of genome is different use win size 100
    try:
        iWinsize = max([int(round(flGenLen/diffs)), 100])
    except ZeroDivisionError:
        return (diffs, d)

    for pos in oBT:
        ref_base = oBT[pos].get("reference")
        s1_base = oBT[pos].get(sample_1, ref_base)
        # consider only differences between valid characters -> pairwise deletion
        if dValChars.get(s1_base.upper(), None) == None:
            continue
        s2_base = oBT[pos].get(sample_2, ref_base)
        # consider only differences between valid characters -> pairwise deletion
        if dValChars.get(s2_base.upper(), None) == None:
            continue
        if s1_base != s2_base:
            iDiffsInWin = 0
            iWinStart = max(0, pos - ((iWinsize / 2) - 1))
            iWinStop = min(int(flGenLen), (pos + (iWinsize / 2) + 1))
            for x in range(iWinStart, iWinStop):
                try:
                    winbase_1 = 'ref' if sample_1 == 'reference' else oBT[x].get(sample_1, 'ref')
                    winbase_2 = 'ref' if sample_2 == 'reference' else oBT[x].get(sample_2, 'ref')
                except KeyError:
                    # x is a point in the alignment where everything is ref
                    continue
                # now x is a point in the alignment where there is at least one SNP but not necessarity in s1 or s2
                if winbase_1 != 'ref':
                    # if not ref check for valid char
                    if dValChars.get(winbase_1.upper(), None) == None:
                        continue
                if winbase_2 != 'ref':
                    # if not ref check for valid char
                    if dValChars.get(winbase_2.upper(), None) == None:
                        continue
                if winbase_1 != winbase_2:
                    iDiffsInWin += 1
            d[pos] = iDiffsInWin

    return (diffs, d)

# --------------------------------------------------------------------------------------------------

def parse_wg_alignment(dArgs, avail_pos, aSampleNames):
    '''
    Parse alignment to data structure
    Parameters
    ----------
    dArgs: dict
        input parameter dictionary as created by get_args()
    avail_pos: dict
        dict of bintrees for each contig
    aSampleNames: list
        list of sample names
    Returns
    -------
    0
    also writes all data to avail_pos
    '''

    """
    avail_pos looks like this:
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

    if os.path.exists(dArgs['alignment_input']) == False:
        logging.error("Input alignment file not found.")
        return 1

    dSeqs = {}

    with open(dArgs['alignment_input'], 'r') as algn:
        # parse the fa file
        sSeq = ""
        sName = ""
        for sLine in algn:
            sLine = sLine.strip()
            if sLine.startswith(">"):
                if len(sSeq) > 0:
                    dSeqs[sName] = sSeq.upper()
                    sSeq = ""
                sName = sLine[1:].split(' ')[0]
                continue
            sSeq = sSeq + sLine
        dSeqs[sName] = sSeq.upper()

    sRefName = 'reference'
    if dSeqs.has_key(sRefName) == False:
        if dArgs['refgenomename'] is not None and dSeqs.has_key(dArgs['refgenomename']) == True:
            sRefName = dArgs['refgenomename']
        else:
            logging.error("No seq named 'reference' in your alignment AND alternative name not found or not given either.")
            return 1
    else:
        logging.info("Your reference appears to be named 'reference'")

    for sn in dSeqs.keys():
        aSampleNames.append(sn)

    all_seq_len = [len(x) for x in dSeqs.values()]
    try:
        assert len(all_seq_len) == all_seq_len.count(all_seq_len[0])
    except AssertionError:
        logging.error("Not all seqs in your alignment have the same length")
        return 1

    avail_pos['alignment_contig'] = FastRBTree()

    iSeqLen = all_seq_len[0]

    for i in range(0, iSeqLen):
        d = {}
        d['reference'] = dSeqs[sRefName][i]
        for sn in aSampleNames:
            if sn != sRefName and dSeqs[sn][i] != d['reference']:
                d[sn] = dSeqs[sn][i]
        # are the nt at this pos all the same?
        if len(d.values()) == d.values().count(d.values()[0]):
            continue
        oBS = BaseStats()
        for sn in aSampleNames:
            oBS.update(d, sn, '_')
        d['stats'] = oBS
        avail_pos['alignment_contig'][i+1] = d

    return 0

# --------------------------------------------------------------------------------------------------

def parse_vcf_files(dArgs, avail_pos, aSampleNames):
    '''
    Parse vcf files to data structure
    Parameters
    ----------
    dArgs: dict
        input parameter dictionary as created by get_args()
    avail_pos: dict
        dict of bintrees for each contig
    aSampleNames: list
        list of sample names
    Returns
    -------
    0
    also writes all data to avail_pos
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
                                                            "stats": BaseStats()})

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
                        logging.info("POS %s passed filters but has multiple alleles REF: %s, ALT: %s. Inserting N", record.POS, str(record.REF), str(record.ALT))
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
    Calculates the distance matrix, optionally removes recombination
    from it and optionally normalises it

    Parameters
    ----------
    aSampleNames: list
        list of sample names
    avail_pos: dict
        infomatin on all available positions
        {'gi|194097589|ref|NC_011035.1|':
            FastRBTree({2329: {'stats': <vcf2distancematrix.BaseStats object at 0x40fb590>,
                               'reference': 'A',
                               '211700_H15498026501': 'C',
                               '211701_H15510030401': 'C',
                               '211702_H15522021601': 'C'},
                        3837: {'211700_H15498026501': 'G',
                               'stats': <vcf2distancematrix.BaseStats object at 0x40fbf90>,
                               '211701_H15510030401': 'G',
                               'reference': 'T',
                               '211702_H15522021601': 'G'},
                        4140: {'211700_H15498026501': 'A',
                               'stats': <vcf2distancematrix.BaseStats object at 0x40fb790>,
                               '211701_H15510030401': 'A',
                               'reference': 'G',
                               '211702_H15522021601': 'A'}})}
    dArgs: dict
        input parameter dictionary as created by get_args()

    Returns
    -------
    call to get_sample_pair_densities with parameters unpacked
    """

    dDen = None
    dRemovals = None
    if dArgs['remove_recombination'] == True:
        if HAVE_SCIPY == False:
            logging.error("Cannot import scipy requied for recombination removal.")
            return None
        dDen = precompute_snp_densities(avail_pos, aSampleNames, dArgs)
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
            else:  # j > i
                pass

    flGenLen = 0.0
    if dDen != None:
        (_, flGenLen) = get_ref_freqs(dArgs['refgenome'], len_only=True)

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

                        # Recombination removal happens here
                        if dDen != None and s1_base != s2_base:

                            iDiffsInWin = dDen[sContig][sample_1][sample_2][pos]
                            iTotalDiffs = dDen['diffs'][sample_1][sample_2]
                            flNofWins = flGenLen / max([int(round(flGenLen/iTotalDiffs)), 100])
                            p_hitting_window = 1.0 / flNofWins
                            p_ok = 1.0

                            # only do binomial test if there are 'too many differences'
                            #  => do not exclude diffs because there are 'not enough'
                            if iDiffsInWin > 1 and iDiffsInWin > iTotalDiffs / float(flNofWins):
                                # binomial test:
                                # what is the probability that I have x successes (i.e. diffs in the current window)
                                # given n trials (i.e. the total number of differences between the two samples)
                                # and given that the probabilty of success is 1/total_num_windows
                                p_ok = binom_test(iDiffsInWin, iTotalDiffs, p_hitting_window)

                            corr_p_thresh = (0.05 / iTotalDiffs)
                            aStats.append("%s\t%s\t%i\t%i\t%i\t%e\t%e\n" % (sample_1,
                                                                            sample_2,
                                                                            pos,
                                                                            iDiffsInWin,
                                                                            iTotalDiffs,
                                                                            p_ok,
                                                                            corr_p_thresh))

                            # null-hypothesis: probability of hitting it is equal for all windows, i.e.
                            # the diffs between the samples are uniformly distributed
                            # Bonferroni corrected p-value threshold: (0.05 / iTotalDiffs)
                            if p_ok <= corr_p_thresh:
                                # likely to be recombinat site (for a given definition of 'likely' and 'recombinant')
                                # aStats.append("%s\t%s\t%i\t%i\t%i\t%e\n" % (sample_1, sample_2, pos, iDiffsInWin, iTotalDiffs, p_ok))
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
                    else:  # j >= i
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
        with open("%s.all_snps.tsv" % (sOutBase), 'w') as fOut:
            fOut.write("sample_1\tsample_2\tposition\tdiffs_in_win\ttotal_diffs\tp_ok\tthreshold\n")
            for sLine in aStats:
                fOut.write(sLine)

    # 'normalise' distance matrix according to model requested
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
    '''
    Normalise distance matrix according to the Tamura 3-parameter distance model
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 16

    Parameters
    ----------
    d: dict
        distance matrix
    ref: str
        reference genome file name
    names: list
        list of sample names
    Returns
    -------
    d: dict
        normalised matrix
    '''

    (dRefFreq, flGenLen) = get_ref_freqs(ref)
    flGC = dRefFreq['C'] + dRefFreq['G']
    h = 2.0 * flGC * (1.0 - flGC)

    for i, sample_1 in enumerate(names):
        for j, sample_2 in enumerate(names):
            if j < i:
                Q = d[sample_1][sample_2][1] / flGenLen
                p = sum(d[sample_1][sample_2]) / flGenLen
                x1 = -h * math.log(1.0 - p / h - Q)
                x2 = 0.5 * (1.0 - h) * math.log(1.0 - (2.0 * Q))
                d[sample_1][sample_2] = x1 - x2
            elif i == j:
                d[sample_1][sample_2] = 0.0

    return d

# --------------------------------------------------------------------------------------------------

def normalise_tn84(d, ref, names):
    """
    Normalise distance matrix according to the Kimura 2-parameter distance model

    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 13 and 14

    Parameters
    ----------
    d: dict
        distance matrix
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
    ref: str
        reference genome file name
    names: list
        list of sample names
    Returns
    -------
    d: dict
        normalised matrix
    """

    (dRefFreq, flGenLen) = get_ref_freqs(ref)

    sum1 = sum([x * x for x in dRefFreq.values()])

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
                if sum2 == 0.0:  # samples are "identical"
                    d[sample_1][sample_2] = 0.0
                else:
                    b = 0.5 * ((1.0 - sum1) + (p * p / sum2))
                    d[sample_1][sample_2] = -b * math.log(1.0 - (p / b))
            elif i == j:
                d[sample_1][sample_2] = 0.0
            else:  # j > i
                pass

    return d

# --------------------------------------------------------------------------------------------------

def normalise_k80(d, ref, names):
    '''
    Normalise distance matrix according to the Tajima-Nei distance model
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 9

    Parameters
    ----------
    d: dict
        distance matrix
    ref: str
        reference genome file name
    names: list
        list of sample names
    Returns
    -------
    d: dict
        normalised matrix
    '''

    (_, flGenLen) = get_ref_freqs(ref, len_only=True)

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
    '''
    Normalise distance matrix according to the Jukes-Cantor distance model
    see: Nei and Zhang: Evolutionary Distance: Estimation,
         ENCYCLOPEDIA OF LIFE SCIENCES 2005,
         doi: 10.1038/npg.els.0005108
         http://www.umich.edu/~zhanglab/publications/2003/a0005108.pdf, equation 7

    Parameters
    ----------
    d: dict
        distance matrix
    ref: str
        reference genome file name
    names: list
        list of sample names
    Returns
    -------
    d: dict
        normalised matrix
    '''

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
                d[sample_1][sample_2] = (-3.0 / 4.0) * math.log(1.0 - ((4.0 / 3.0) * p))

    return d

# --------------------------------------------------------------------------------------------------

def get_difference_value(s1_base, s2_base, sSubs):
    '''
    Get difference value for a given set of bases.

    Parameters
    ----------
    s1_base: str
        a charcater
    s2_base: str
        a charcater
    sSubs: str
        distance model

    Returns
    -------
    difference: float or list
        depending on the distance model either a float 1.0 or 0.0
        of a list of two floats
    '''
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
    Sum up total number of differences for a dict like the one in the input

    Parameters
    ----------
    d: dict
        {'A': {'A': 1152.0, 'C': 114.0, 'G': 545.0, 'T': 35.0},
         'C': {'A': 0.0, 'C': 1233.0, 'G': 108.0, 'T': 467.0},
         'G': {'A': 0.0, 'C': 0.0, 'G': 1283.0, 'T': 100.0},
         'T': {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 1177.0}}

    Returns
    -------
    t: float
        the sum of all differences(1369.0 in above example case)
    """

    t = 0.0
    for k, v in d.items():
        t += sum([y if x != k else 0.0 for x, y in v.items()])
    return t

# --------------------------------------------------------------------------------------------------

def get_ref_freqs(ref, len_only=False):
    """
    Get the length of the reference genome and optionally the nucleotide frequencies in it

    Parameters
    ----------
    ref: str
        reference genome filename
    len_only: boolean
        get genome lengths only [default FALSE, also get nucleotide frequencies]

    Returns
    -------
    (dRefFreq, flGenLen): tuple
        dRefFreq: dict
            dRefFreq = {'A': 0.25, 'C': 0.24, 'G': 0.26, 'T': 0.25}
        flGenLen: float
            genome length
    """

    # get frequency of a, c, g, and t in ref
    dRefFreq = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
    flGenLen = 0.0
    with open(ref, 'r') as fRef:
        for sLine in fRef:
            if sLine.startswith(">"):
                continue
            else:
                sLine = sLine.strip()
                flGenLen += len(sLine)
                if len_only == True:
                    continue
                sLine = sLine.upper()
                for n in sLine:
                    try:
                        dRefFreq[n] += 1.0
                    except KeyError:
                        pass
    if len_only == False:
        for (i, j) in dRefFreq.items():
            dRefFreq[i] = j / flGenLen

    logging.info("genome length: %s" % flGenLen)
    return (dRefFreq, flGenLen)

# --------------------------------------------------------------------------------------------------

def calculate_memory_for_sort():
    """Calculate available memory for ``samtools sort`` function.
    If there is enough memory, no temp files are created. **Enough**
    is defined as at least 1G per CPU.

    Returns
    -------
    sort_memory: str or None
        String to use directly with *-m* option in sort, or None.
    """

    avail_memory = None
    if HAVE_PSUTIL == True:
        mem = virtual_memory()
        avail_memory = mem.total
    else:
        try:
            avail_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
        except ValueError:
            logging.error("If you're in Mac OS you need to have the psutil Python library.")
            raise SystemExit

    avail_cpu = multiprocessing.cpu_count()

    sort_memory = avail_memory / avail_cpu / 1024 ** 2

    # samtools default from documentation is 768M, to be conservative
    #    only use -m option when there is more than 1G per CPU.
    if sort_memory < 1024:
        sort_memory = None

    else:
        sort_memory = "%sG" % (int(floor(sort_memory / 1024)))

    return sort_memory

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    print calculate_memory_for_sort()
