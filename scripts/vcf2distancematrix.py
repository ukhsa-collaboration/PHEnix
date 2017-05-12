#!/usr/bin/env python
'''
Provides the vcf to diastance matrix funtionality plus the option to remove
differences between samples unlikely to be cause by background drift (recombination removal)
'''

import argparse
import glob
import logging
import os
import sys

from Bio import Phylo
from Bio.Phylo import TreeConstruction

from phe.utils import parse_vcf_files, get_dist_mat, parse_wg_alignment


# --------------------------------------------------------------------------------------------------
def get_desc():
    """
    Get description for help text

    Parameters
    ----------
    no inputs

    Returns
    -------
    a string with contains the description
    """

    return """Combine multiple VCFs into a distance matrix.
              Distance measures according to five different models are available:
              * Number of differences\n
              * Jukes-Cantor distance (jc69)
              * Tajima-Nei distance (k80)
              * Kimura 2-parameter distance (tn84)
              * Tamura 3-parameter distance (t93)
           """

# --------------------------------------------------------------------------------------------------

def get_args():
    """
    Parse arguments

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
                       help="Path to the directory with .vcf files. Input option 1.")

    group.add_argument("--input",
                       "-i",
                       type=str,
                       nargs='+',
                       help="List of VCF files to process. Input option 2.")

    group.add_argument("--alignment-input",
                       "-a",
                       metavar="MULTI FASTA FILE",
                       help="Multi fasta file with whole genome input alignment. Input option 3.")

    parser.add_argument("--out",
                        "-o",
                        required=True,
                        help="Path to the maxtrix output file in given format. \
                              [REQUIRED. default format is tab separated. \
                              use --format to change format]")

    parser.add_argument('--deletion',
                        metavar="STRING",
                        default='pairwise',
                        choices=['pairwise', 'complete'],
                        dest="deletion",
                        help="Method of recombination filtering. Either 'pairwise' \
                              or 'complete' ['pairwise']")

    parser.add_argument('--substitution',
                        metavar="STRING",
                        default='number_of_differences',
                        choices=['number_of_differences', 'jc69', 'k80', 'tn84', 't93'],
                        dest="substitution",
                        help="Substituition model. Either 'number_of_differences', \
                              'jc69', 'k80', 'tn84' or 't93' ['number_of_differences']")

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
                        help="Reference genome used for SNP calling [Required for recombination removal].")

    parser.add_argument("--refgenomename",
                        "-n",
                        type=str,
                        metavar="STRING",
                        dest="refgenomename",
                        default=None,
                        help="Name of reference genome in input alignment [Required if input option 3 is used and reference is not named 'reference'].")

    parser.add_argument("--threshold",
                        "-k",
                        type=float,
                        metavar="FLOAT",
                        dest="k",
                        default=1.0,
                        help="Density tyhreshold above mean density for relevant pair. [1.0].")

    parser.add_argument("--threads",
                        type=int,
                        metavar="INT",
                        dest="threads",
                        default=1,
                        help="Number of threads to use. [1].")

    parser.add_argument("--format",
                        type=str,
                        metavar="STRING",
                        dest="format",
                        choices=['tsv', 'csv', 'mega'],
                        default='tsv',
                        help="Change format for output file. Available \
                              options csv, tsv and mega. [tsv]")

    parser.add_argument("--tree",
                        "-t",
                        type=str,
                        metavar="FILE",
                        dest="tree",
                        default=None,
                        help="Make an NJ tree and write it to the given file in newick \
                              format. [Default: Don't make tree, only matrix]")

    parser.add_argument('--with-stats',
                        action="store_true",
                        default=False,
                        help="Write additional files with information on \
                              removed recombinant SNPs. [don't]")

    return parser

# --------------------------------------------------------------------------------------------------

def main(dArgs):
    '''
    Main funtion

    Parameters
    ----------
    dArgs: dict
        input parameter dictionary as created by get_args()

    Returns
    -------
    returns 0
    also creates all logs and result files
    '''

    if dArgs['substitution'] in ['jc69', 'k80', 'tn84', 't93'] and dArgs['refgenome'] == None:
        sys.stderr.write("Error: Please specify reference genome for substitution model.\n")
        return 1

    if dArgs['remove_recombination'] == True:
        try:
            from scipy.stats import binom_test
        except ImportError:
            sys.stderr.write("ERROR: cannot import scipy which is required for recombinatin removal. Exiting!\n")
            sys.stderr.write("Hint: You can still get a distance matrix without recombination removed.\n")
            return 1

    aSampleNames = []
    avail_pos = {}

    if dArgs['alignment_input'] is None:

        if dArgs['directory'] is not None and dArgs['input'] is None:
            dArgs['input'] = glob.glob(os.path.join(dArgs['directory'], "*.vcf"))

        if len(dArgs['input']) <= 0:
            sys.stderr.write("Error: No VCFs found.\n")
            return 1
        else:
            logging.info("%i VCFs found", len(dArgs['input']))
            # parse vcf files into avail_pos structure
            parse_vcf_files(dArgs, avail_pos, aSampleNames)
    else:
        logging.info("Submitted input as multi fasta aligment")
        parse_wg_alignment(dArgs, avail_pos, aSampleNames)

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

    number_of_sites = sum([len(x) for _, x in avail_pos.items()])
    logging.info("%i total variant positions found" % (number_of_sites))

    if dArgs['deletion'] == 'complete':
        for _, oBT in avail_pos.items():
            to_del = []
            for iPos in oBT:
                if oBT[iPos]['stats'].N > 0 or oBT[iPos]['stats'].gap > 0:
                    to_del.append(iPos)
            for r in to_del:
                oBT.remove(r)
        logging.info("%i total variant positions left after complete removal",
                     sum([len(x) for _, x in avail_pos.items()]))
    else:  # deletion is pairwise, which is implicit during matrix creation
        pass

    dist_mat = {}
    dist_mat = get_dist_mat(aSampleNames, avail_pos, dArgs)

    if dist_mat == None:
        logging.info("An error occured during matrix creation.")
        return 1

    if dArgs['format'] == 'mega':
        write_mega_file(dArgs, aSampleNames, dist_mat, number_of_sites)
    else:
        sep = '\t' if dArgs['format'] == 'tsv' else ','
        sform = "%s%e"
        if dArgs['substitution'] == 'number_of_differences':
            sform = "%s%i"
        with open(dArgs['out'], "wb") as fp:
            for i, sample_1 in enumerate(aSampleNames):
                row = sample_1
                for j, sample_2 in enumerate(aSampleNames):
                    if j < i:
                        dist = dist_mat[sample_1][sample_2]
                        row += sform % (sep, dist)
                fp.write("%s\n" % row)

    if dArgs['tree'] != None:
        make_nj_tree(dist_mat, dArgs, aSampleNames)

    logging.info("Done.")

    return 0

# --------------------------------------------------------------------------------------------------

def write_mega_file(dArgs, aSampleNames, dist_mat, number_of_sites=0):
    '''
    Writes a mega output file of a distance matrix.

    Parameters
    ----------
    dArgs: dict
        input argument dictionary
    aSampleNames: list
        list of sample names
    dist_mat: dict
        distance matrix as a dict of dicts
        distance_a_to_b = dist_mat[a][b]
    number_of_sites: int
        number of sites to put in mega header (default 0)

    Returns
    -------
    returns 0
    also writes files to dArgs['out']
    '''

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
            fp.write("[%s] #%s\n" % (str(i + 1).rjust(spacing1), name))

        fp.write('\n')
        fp.write('[              %s ]\n' \
                 % ('            '.join([str(x) for x in range(1, len(aSampleNames) + 1)])))
        for i, sample_1 in enumerate(aSampleNames):
            row = "[%s] " % (str(i + 1).rjust(spacing1))
            for j, sample_2 in enumerate(aSampleNames):
                if j < i:
                    dist = dist_mat[sample_1][sample_2]
                    row += " %12.10f" % (dist)
            fp.write("%s\n" % row)

    return 0

# end of main --------------------------------------------------------------------------------------

def make_nj_tree(dist_mat, dArgs, aSampleNames):
    '''
    Uses Biopython.Phylo to make a neighbour joining tree from a distance matrix

    Parameters
    ----------
    dist_mat: dict
        distance matrix as a dict of dicts
        distance_a_to_b = dist_mat[a][b]
    dArgs: dict
        input argument dictionary
    aSampleNames: list
        list of sample names

    Returns
    -------
    returns 0
    also writes tree file to to dArgs['tree'] in newick format
    '''

    aSimpleMatrix = []
    for i, sample_1 in enumerate(aSampleNames):
        mat_line = []
        for j, sample_2 in enumerate(aSampleNames):
            if j < i:
                mat_line.append(dist_mat[sample_1][sample_2])
            elif j == i:
                mat_line.append(0)
            else:
                pass
        aSimpleMatrix.append(mat_line)

    oDistMat = TreeConstruction._DistanceMatrix(aSampleNames, aSimpleMatrix)
    constructor = TreeConstruction.DistanceTreeConstructor()
    oTree = constructor.nj(oDistMat)
    Phylo.write(oTree, dArgs['tree'], 'newick')
    logging.info("Tree file written.")

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    exit(main(vars(get_args().parse_args())))
