#!/usr/bin/env python
'''
Create a json representation of the variant and ignored positions in a VCF file and write to a file

:Date: 28 November 2016
:Author: Anthony Underwood
'''

import sys
import argparse
import vcf
import gzip
import json
import logging

from phe.utils import parse_vcf_files

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

    return """Converts the postions of variants and ignored/missing positions in either a 'raw' or filtered VCF
              file to a json string and writes it to a file.
              The json contains 6 arrays for each chromosome in the VCF file:
              g_positions, a_positions, t_positions, c_positions, gap_positions, n_positions
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

    parser.add_argument("--input",
                       "-i",
                       required=True,
                       type=str,
                       help="path to a VCF file")

    parser.add_argument("--output_file_prefix",
                        "-o",
                        default='sample_name',
                        type=str,
                        help="Path to the json output file (without file extension). Default: sample_name")

    parser.add_argument('--nozip',
                        "-z",
                        action="store_true",
                        help="Do not gzip json when writing file. (default: Yes, gzip it.)")

    parser.add_argument('--vcf_is_filtered',
                        "-f",
                        action="store_true",
                        help="Required: Confirm that the input vcf is filtered. It is strongly" + \
                             "recommended to filter the file with Phenix using the same" + \
                             "parameters that are used throughout the database this jason file is meant for.")

    parser.add_argument('--summary_info',
                        "-s",
                        action="store_true",
                        help="Print summary of the json string")

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

    if dArgs['vcf_is_filtered'] == False:
        logging.error("Please use the --vcf_is_filtered option to confirm that your vcf has been filtered!")
        logging.error("See the helptext for this option for further info.\nExiting.")
        return 1


    data = {}
    data['annotations'] = {}
    annotations = []

    with open(dArgs['input']) as fVcf:
        sLine = fVcf.readline().strip()
        while sLine.startswith("##") == True:
            if not (sLine.startswith("##INFO") or sLine.startswith("##FORMAT") or sLine.startswith("##FILTER")):
                annotations.append(sLine)
            sLine = fVcf.readline().strip()

    for annot in annotations:
        [annot_name, annot_value] = annot.split("=", 1)
        annot_name = annot_name.replace("#", "")
        annot_value = annot_value.replace('>', '').replace('<', '')
        data['annotations'][annot_name] = annot_value

    #parse_vcf_files expects dArgs to look a little different
    dArgs["include"] = None
    dArgs["exclude"] = None
    dArgs["input"] = [dArgs["input"]]

    aSampleNames = []
    avail_pos = {}

    # parse vcf files into avail_pos structure
    parse_vcf_files(dArgs, avail_pos, aSampleNames)

    sname = filter(lambda x: x != 'reference', aSampleNames)[0]

    data['positions'] = {}

    for chrom, rbtree in avail_pos.iteritems():
        data['positions'][chrom] = {"G": [], "A": [], "T": [], "C": [], "N": [], "-": []}
        for pos, pos_info in rbtree.items():
            sample_base = pos_info[sname]
            if sample_base != pos_info['reference']:
                data['positions'][chrom][sample_base].append(pos)

    # print summary information if requested
    if dArgs["summary_info"]:
        logging.info("%s chromosomes were found in the VCF file" % len(data['positions']))
        for chromosome in data['positions']:
            logging.info("")
            logging.info("%s: " % chromosome)
            for position_type in ["G", "A", "T", "C", "N", "-"]:
                logging.info("\t%s: %s" % (position_type, len(data['positions'][chromosome][position_type])))
            logging.info("-----------------")

    out_file_prefix = dArgs["output_file_prefix"] if dArgs["output_file_prefix"] != 'sample_name' else sname

    # create json string
    json_string = json.dumps(data)
    if dArgs["nozip"] == False:
        with gzip.open("%s.json.gz" % (out_file_prefix), mode='wb') as gzip_obj:
            gzip_obj.write(json_string)
    else:
        with open("%s.json" % (out_file_prefix), mode='w') as file_obj:
            file_obj.write(json_string)

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    exit(main(vars(get_args().parse_args())))
