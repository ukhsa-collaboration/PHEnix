#!/usr/bin/env python
'''

'''
import sys
import argparse
import glob
import logging
import os
import pprint

from Bio import Phylo
from Bio.Phylo import TreeConstruction

from phe.utils import base_stats, parse_vcf_files
from phe.utils import is_uncallable, get_dist_mat

# --------------------------------------------------------------------------------------------------

def get_desc():
    return """Combine multiple VCFs into a distance matrix.
              Distance measures according to five different models are available:
              * Number of differences\n
              * Page Jukes-Cantor distance (jc69)
              * Tajima-Nei distance (k80)
              * Kimura 2-parameter distance (tn84)
              * Tamura 3-parameter distance (t93)
           """

# --------------------------------------------------------------------------------------------------

def get_args():
    """
    Parge arguments
    Parameters
    ----------
    no inputs
    Returns
    -------
    oArgs: obj
        arguments object
    """

    parser = argparse.ArgumentParser(description=get_desc())

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("--directory",
                       "-d",
                       help="Path to the directory with .vcf files.")

    group.add_argument("--input",
                       "-i",
                       type=str,
                       nargs='+',
                       help="List of VCF files to process.")

    parser.add_argument("--out",
                        "-o",
                        required=True,
                        help="Path to the maxtrix output file in given format. [REQUIRED. default format is tab separated. use --format to change format]")

    parser.add_argument('--deletion',
                        metavar="STRING",
                        default='pairwise',
                        choices=['pairwise', 'complete'],
                        dest="deletion",
                        help="Method of recombination filtering. Either 'pairwise' or 'complete' ['pairwise']")

    parser.add_argument('--substitution',
                        metavar="STRING",
                        default='number_of_differences',
                        choices=['number_of_differences', 'jc69', 'k80', 'tn84', 't93'],
                        dest="substitution",
                        help="Substituition model. Either 'number_of_differences', 'jc69', 'k80', 'tn84' or 't93' ['number_of_differences']")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--include",
                       metavar="BED FILE",
                       dest='include',
                       default=None,
                       help="Only include positions in BED file in the FASTA")

    group.add_argument("--exclude",
                       metavar="BED FILE",
                       dest='exclude',
                       default=None,
                       help="Exclude any positions specified in the BED file.")

    parser.add_argument('--remove-recombination',
                        action="store_true",
                        default=False,
                        help="Attempt to remove recombination from distance matrix. [don't]")

    parser.add_argument("--refgenome",
                        "-g",
                        type=str,
                        metavar="FASTA FILE",
                        dest="refgenome",
                        default=None,
                        help="Reference genome used for SNP calling [Required for 'jc69', 'k80', 'tn84' and 't93' substitution, else ignored].")

    parser.add_argument("--window-size",
                        "-w",
                        type=int,
                        metavar="INTEGER",
                        dest="winsize",
                        default=200,
                        help="Window size of SNP density calculation. [200]")

    parser.add_argument("--max-snps",
                        "-m",
                        type=int,
                        metavar="INTEGER",
                        dest="maxsnps",
                        default=5,
                        help="Maximum number of SNPs in window with no SNPs in other sample [5].")

    parser.add_argument("--threshold",
                        "-r",
                        type=float,
                        metavar="FLOAT",
                        dest="threshold",
                        default=2.0,
                        help="Maxium allowable density ratio in two windows.[2.0 i.e. double].")

    parser.add_argument("--format",
                        type=str,
                        metavar="STRING",
                        dest="format",
                        choices=['tsv', 'csv', 'mega'],
                        default='tsv',
                        help="Change format for output file. Available options csv and tsv.")

    return parser

# --------------------------------------------------------------------------------------------------

def main(dArgs):
    '''
    Main funtion
    Parameters
    ----------
    no inputs
    Returns
    -------
    0
    Creates all logs and result files
    '''

    if dArgs['substitution'] in ['jc69', 'k80', 'tn84', 't93'] and dArgs['refgenome'] == None:
        sys.stderr.write("Error: Please specify reference genome for substitution model.\n")
        return 1

    if dArgs['directory'] is not None and dArgs['input'] is None:
        dArgs['input'] = glob.glob(os.path.join(dArgs['directory'], "*.vcf"))

    if len(dArgs['input']) <= 0:
        sys.stderr.write("Error: No VCFs found.\n")
        return 1
    else:
        logging.info("%i VCFs found" %(len(dArgs['input'])))

    aSampleNames = []
    avail_pos = {}

    # parse vcf files into avail_pos structure
    parse_vcf_files(dArgs, avail_pos, aSampleNames)

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

    number_of_sites = sum([len(x) for _, x in avail_pos.items()])
    logging.info("%i total variant positions found" % (number_of_sites))

    if dArgs['deletion'] == 'complete':
        for contig, oBT in avail_pos.items():
            for iPos in oBT:
                if oBT[iPos]['stats'].N > 0 or oBT[iPos]['stats'].gap > 0:
                    oBT.remove(iPos)
        logging.info("%i total variant positions left after complete removal" % (sum([len(x) for _, x in avail_pos.items()])))
    else: # deletion is pairwise, which is implicit during matrix creation
        pass

    dist_mat = {}
    dist_mat = get_dist_mat(aSampleNames, avail_pos, dArgs)

    pprint.pprint(dist_mat)


    if dArgs['format'] == 'mega':
        write_mega_file(dArgs, aSampleNames, dist_mat, number_of_sites)
    else:
        sep = '\t' if dArgs['format'] == 'tsv' else ','
        # aSimpleMatrix = []
        with open(dArgs['out'], "wb") as fp:
            # fp.write("%s%s\n" % (sep, sep.join(aSampleNames)))
            for i, sample_1 in enumerate(aSampleNames):
                row = sample_1
                # mat_line = []
                for j, sample_2 in enumerate(aSampleNames):
                    if j < i:
                        dist = dist_mat[sample_2][sample_1]
                        # mat_line.append(dist)
                    # elif j == i:
                        # mat_line.append(0)
                    #     dist = dist_mat[sample_1][sample_2]
                    # else:
                    #    dist = dist_mat[sample_1][sample_2]
                        row += "%s%e" % (sep, dist)
                fp.write("%s\n" % row)
                # aSimpleMatrix.append(mat_line)

    # oDistMat = TreeConstruction._DistanceMatrix(aSampleNames, aSimpleMatrix)
    # constructor = TreeConstruction.DistanceTreeConstructor()
    # oTree = constructor.nj(oDistMat)
    # Phylo.write(oTree, oArgs.tree,'newick')
    # if oArgs.ascii == True:
    #     Phylo.draw_ascii(oTree)

    return 0

# --------------------------------------------------------------------------------------------------

def write_mega_file(dArgs, aSampleNames, dist_mat, number_of_sites=0):

    header = """#mega
!Title: fasta file;
!Format DataType=Distance DataFormat=LowerLeft NTaxa=<number_of_taxa>;
!Description
  Analysis
    Analysis ---------------------- Distance Estimation
    Scope ------------------------- Pairs of taxa
  Estimate Variance
    Variance Estimation Method ---- None
  Substitution Model
    Substitutions Type ------------ Nucleotide
    Model/Method ------------------ <model>
    Substitutions to Include ------ d: Transitions + Transversions
  Rates and Patterns
    Rates among Sites ------------- Uniform rates
    Pattern among Lineages -------- Same (Homogeneous)
  Data Subset to Use
    Gaps/Missing Data Treatment --- <deletion> deletion
  No. of Sites : <number_of_sites>
  d : Estimate
;"""

    header = header.replace('<number_of_taxa>', str(len(aSampleNames)))
    header = header.replace('<model>', dArgs['substitution'])
    header = header.replace('<deletion>', dArgs['deletion'])
    header = header.replace('<number_of_sites>', str(number_of_sites))

    spacing1 = len(str(len(aSampleNames)))
    with open(dArgs['out'], "wb") as fp:
        fp.write('%s\n\n' % (header))
        for i, name in enumerate(aSampleNames):
            fp.write("[%s] #%s\n" % (str(i+1).rjust(spacing1), name))

        fp.write('\n')
        fp.write('[              %s ]\n' % ('            '.join([str(x) for x in range(1, len(aSampleNames) + 1)])))
        for i, sample_1 in enumerate(aSampleNames):
            row = "[%s] " % (str(i+1).rjust(spacing1))
            for j, sample_2 in enumerate(aSampleNames):
                if j < i:
                    dist = dist_mat[sample_2][sample_1]
                    row += " %12.10f" % (dist)
            fp.write("%s\n" % row)




    return 0


# end of main --------------------------------------------------------------------------------------






















def get_dist_mat_average(aSampleNames, avail_pos, dValChars, oArgs):

    avg_den_mat = {}
    dist_mat = {}
    for i, sample_1 in enumerate(aSampleNames):
        avg_den_mat[sample_1] = {}
        dist_mat[sample_1] = {}
        for j, sample_2 in enumerate(aSampleNames):
            if j < i:
                continue
            avg_den_mat[sample_1][sample_2] = 0
            dist_mat[sample_1][sample_2] = 0

    for i, sample_1 in enumerate(aSampleNames):
        for j, sample_2 in enumerate(aSampleNames):
            if j <= i:
                continue
            iNofSNPs = 0
            for sContig in avail_pos.keys():
                for pos in avail_pos[sContig]:
                    ref_base = avail_pos[sContig][pos].get("reference").upper()
                    s1_base = avail_pos[sContig][pos].get(sample_1, ref_base).upper()
                    s2_base = avail_pos[sContig][pos].get(sample_2, ref_base).upper()
                    if dValChars.get(s1_base, None) == None:
                        continue
                    if dValChars.get(s2_base, None) == None:
                        continue
                    if s1_base != s2_base:
                        iNofSNPs += 1
            avg_den_mat[sample_1][sample_2] = iNofSNPs/float(oArgs.genomesize)

    for sContig in avail_pos.keys():
        for pos in avail_pos[sContig]:
            ref_base = avail_pos[sContig][pos].get("reference")
            # Position has been seen or no reference available.
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
                    if s1_base != s2_base:
                        # dist_mat[sample_1][sample_2] += 1
                        if sample_1 == 'reference' or sample_2 == 'reference':
                            dist_mat[sample_1][sample_2] += 1
                        else:
                            if is_this_a_recombinant_site_average(avail_pos,
                                                                  sContig,
                                                                  pos,
                                                                  sample_1,
                                                                  sample_2,
                                                                  avg_den_mat[sample_1][sample_2],
                                                                  oArgs.winsize,
                                                                  oArgs.threshold,
                                                                  dValChars) == False:
                                dist_mat[sample_1][sample_2] += 1

    return dist_mat

# --------------------------------------------------------------------------------------------------

def get_dist_mat_pairwise(aSampleNames, avail_pos, dValChars, oArgs):

    dist_mat = {}
    for i, sample_1 in enumerate(aSampleNames):
        dist_mat[sample_1] = {}
        for j, sample_2 in enumerate(aSampleNames):
            if j < i:
                continue
            dist_mat[sample_1][sample_2] = 0

    for sContig in avail_pos.keys():
        for pos in avail_pos[sContig]:
            ref_base = avail_pos[sContig][pos].get("reference")
            # Position has been seen or no reference available.
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
                    if s1_base != s2_base:
                        # dist_mat[sample_1][sample_2] += 1
                        if sample_1 == 'reference' or sample_2 == 'reference':
                            dist_mat[sample_1][sample_2] += 1
                        else:
                            if is_this_a_recombinant_site(avail_pos,
                                                          sContig,
                                                          pos,
                                                          sample_1,
                                                          sample_2,
                                                          oArgs.winsize,
                                                          oArgs.maxsnps,
                                                          oArgs.threshold) == False:
                                dist_mat[sample_1][sample_2] += 1

    return dist_mat

# --------------------------------------------------------------------------------------------------

def is_this_a_recombinant_site_average(avail_pos,
                                       sContig,
                                       pos,
                                       sample_1,
                                       sample_2,
                                       avg_den,
                                       winsize,
                                       threshold,
                                       dValChars):

    bIsIt = False

    den = 0.0
    half_window = int((winsize/2.0)+0.5)

    start = pos - half_window
    start = 0 if start < 0 else start

    snp_cnt = 0.0
    pos_info = None
    for i in range(start, pos + half_window):

        try:
            pos_info = avail_pos[sContig][i]
        except KeyError:
            continue

        ref_base = pos_info.get("reference").upper()
        s1_base = pos_info.get(sample_1, ref_base).upper()
        s2_base = pos_info.get(sample_2, ref_base).upper()
        if dValChars.get(s1_base, None) == None:
            continue
        if dValChars.get(s2_base, None) == None:
            continue
        if s1_base != s2_base:
            snp_cnt += 1.0

    den = snp_cnt/(2.0*half_window)

    if den > avg_den * threshold:
        bIsIt = True

    # print snp_cnt, den, avg_den, winsize, threshold, bIsIt

    return bIsIt

# --------------------------------------------------------------------------------------------------

def is_this_a_recombinant_site(avail_pos, contig, pos, sample_1, sample_2, window_size, max_snps, threshold):

    bIsIt = False

    (density_1, snp_cnt_1) = get_local_snp_density(avail_pos, contig, pos, sample_1, window_size)
    (density_2, snp_cnt_2) = get_local_snp_density(avail_pos, contig, pos, sample_2, window_size)

    assert snp_cnt_1 > 0.0 or snp_cnt_2 > 0.0

    den_ratio = 0.0

    if snp_cnt_2 == 0.0 or snp_cnt_1 == 0.0:
        if max(snp_cnt_1, snp_cnt_2) >= float(max_snps):
            bIsIt = True
        else:
            bIsIt = False
    else:
        den_ratio = density_1 / density_2

        if den_ratio >= (1/float(threshold)) and den_ratio <= float(threshold):
            bIsIt = False
        else:
            bIsIt = True

    return bIsIt

# --------------------------------------------------------------------------------------------------

def get_local_snp_density(avail_pos, contig, pos, sample, window_size):

    den = 0.0
    half_window = int((window_size/2.0)+0.5)

    start = pos - half_window
    start = 0 if start < 0 else start

    snp_cnt = 0.0
    pos_info = None
    for i in range(start, pos + half_window):

        try:
            pos_info = avail_pos[contig][i]
        except KeyError:
            continue

        snp_base = pos_info.get(sample, None)
        if snp_base != None and snp_base.upper() in ["A", "C", "G", "T"]:
            snp_cnt += 1.0

    den = snp_cnt/(2.0*half_window)

    return (den, snp_cnt)

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    exit(main())
