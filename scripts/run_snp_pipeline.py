from argparse import RawTextHelpFormatter
import argparse
from collections import OrderedDict
import glob
import logging
import os
import sys

import yaml

from phe.annotations import available_annotators, make_annotators
from phe.mapping.mapping_factory import factory as map_fac, available_mappers
from phe.variant import VariantSet
from phe.variant.variant_factory import factory as variant_fac, \
    available_callers
from phe.variant_filters import available_filters, str_to_filters, make_filters


def get_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "VERSION")
    version = "N/A"

    if os.path.exists(version_file):
        try:
            with open(version_file) as fp:
                version = fp.next().strip()
                version += "-static"
        except IOError:
            pass
    return version

def pipeline(workflow, input_dir):
    '''Setup config for pipeline run.

    Parameters:
    -----------
    workflow: str
        Name of the workflow running.
    input_dir: str
        Path to the input directory where the fastQ files are kept.

    Returns:
    --------
    dict:
        Dictionary with r1, r2, reference, sample_name, config and outdir keys.
    '''
    config = {}
    for fastq in glob.glob(os.path.join(input_dir, "*processed*fastq.gz")):
        if fastq.endswith("R1.fastq.gz"):
            config["r1"] = os.path.abspath(fastq)

            config["sample_name"] = ".".join(os.path.basename(fastq).split(".")[:-6])
        elif fastq.endswith("R2.fastq.gz"):
            config["r2"] = os.path.abspath(fastq)

    output_dir = os.path.join(input_dir, "snp_pipeline")

    if not os.path.exists(output_dir):
        # os.makedirs(output_dir)
        logging.critical("Output directory %s doesn't exists, Exiting.", output_dir)
    else:
        config["outdir"] = output_dir

    config_files = glob.glob(os.path.join(output_dir, "*.yml"))
    if len(config_files) > 1:
        logging.critical("More than 1 YAML file detected in %s. A single config is required.", output_dir)
    elif len(config_files) == 0:
        logging.critical("No YAML configs detected in %s. A config is required.", output_dir)
    else:
        config["config"] = config_files[0]

    reference_file = os.path.join(output_dir, "reference.fasta")

    if not os.path.exists(reference_file):
        logging.critical("No reference found in %s. A reference.fasta is required.", output_dir)
    else:
        config["reference"] = reference_file

    return config

def get_desc():
    return r'''Run the snp pipeline with specified mapper, variant caller and some filters.

Available mappers: %s

Available variant callers: %s

Available filters: %s

Available annotators: %s''' % (available_mappers(), available_callers(), available_filters(), available_annotators())

def get_args():
    args = argparse.ArgumentParser(description=get_desc(), formatter_class=RawTextHelpFormatter)

    args.add_argument("--workflow", "-w")
    args.add_argument("--input", "-i")

    args.add_argument("-r1", help="R1/Forward read in Fastq format.")
    args.add_argument("-r2", help="R2/Reverse read in Fastq format.")
    args.add_argument("--reference", "-r", help="Rerefence to use for mapping.")
    args.add_argument("--sample-name", default="test_sample", help="Name of the sample for mapper to include as read groups.")
    args.add_argument("--outdir", "-o")

    args.add_argument("--config", "-c")

    args.add_argument("--mapper", "-m", default="bwa", help="Available mappers: %s" % available_mappers())
    args.add_argument("--mapper-options", help="Custom maper options (advanced)")
    args.add_argument("--bam")
    args.add_argument("--variant", "-v", default="gatk", help="Available variant callers: %s" % available_callers())
    args.add_argument("--variant-options", help="Custom variant options (advanced)")
    args.add_argument("--vcf")
    args.add_argument("--filters", type=str, help="Filters to be applied to the VCF in key:value pairs, separated by comma (,). Available_filters: %s. Recommendations: GATK: mq_score:30,min_depth:10,ad_ratio:0.9 Mpileup: mq_score:30,min_depth:10,dp4_ratio:0.9" % available_filters())

    args.add_argument("--annotators", nargs="+", help="List of annotators to run before filters. Available: %s" % available_annotators())

    args.add_argument("--keep-temp", action="store_true", help="Keep intermediate files like BAMs and VCFs (default: False).")

    args.add_argument("--json", action="store_true", help="Also write variant positions in filtered vcf as json file (default: False).")
    args.add_argument("--json-info", action="store_true", help="When writing a json file, log some stats to stdout. (default: False).")

    return args

def load_config(args):

    with open(args["config"]) as fp:
        config = yaml.load(fp)

    args["mapper"] = config.get("mapper")
    args["mapper_options"] = config.get("mapper-options")

    args["variant"] = config.get("variant")
    args["variant_options"] = config.get("variant-options")

    args["filters"] = config.get("filters")

    args["annotators"] = config.get("annotators")


def main(args):

    if args.get("version") is None:
        args["version"] = get_version()

    make_aux = False
    if args["workflow"] and args["input"]:
        logging.info("PIPELINE_START")
        workflow_config = pipeline(args["workflow"], args["input"])
        try:
            args["r1"] = workflow_config["r1"]
            args["r2"] = workflow_config["r2"]
            args["outdir"] = workflow_config["outdir"]
            args["config"] = workflow_config["config"]
            args["reference"] = workflow_config["reference"]
            args["sample_name"] = workflow_config["sample_name"]
            make_aux = True
        except KeyError:
            logging.critical("Could not find parameters in %s", args["input"])
            return 5

    logging.info("Initialising data matrix.")

    if args["outdir"] is None:
        sys.stdout.write("Please provide output directory.")
        return -1
    elif not os.path.exists(args["outdir"]):
        os.makedirs(args["outdir"])

    # If config is specified, then load data from that.
    if args["config"]:
        load_config(args)

    mapper = None
    if args["mapper"]:
        mapper = map_fac(mapper=args["mapper"], custom_options=args["mapper_options"])

    variant = None
    if args["variant"]:
        variant = variant_fac(variant=args["variant"], custom_options=args["variant_options"])

    if args["annotators"]:
        args["annotators"] = make_annotators(args["annotators"])

    if args["filters"]:
        try:
            if isinstance(args["filters"], str):
                args["filters"] = str_to_filters(args["filters"])
            elif isinstance(args["filters"], dict):
                args["filters"] = make_filters(args["filters"])
            else:
                logging.warn("Unknown filters specified: %s", args["filters"])
        except Exception:
            logging.error("Failed to recognise and create filters.")
            return 3

    if args["bam"] is not None:
        logging.info("Found BAM file: %s", args["bam"])
        bam_file = args["bam"]
    elif args["vcf"] is None and mapper is not None:
        logging.info("Mapping data file with %s.", args["mapper"])
        bam_file = os.path.join(args["outdir"], "%s.bam" % args["sample_name"])
        success = mapper.make_bam(ref=args["reference"],
                                  R1=args["r1"],
                                  R2=args["r2"],
                                  out_file=bam_file,
                                  sample_name=args["sample_name"],
                                  make_aux=make_aux)

        if not success:
            logging.warn("Could not map reads to the reference. Aborting.")
            return 1
    else:
        bam_file = None

    logging.info("Creating digitised variants with %s.", args["variant"])
    if args["vcf"]:
        vcf_file = args["vcf"]
    elif bam_file is not None:
        vcf_file = os.path.join(args["outdir"], "%s.vcf" % args["sample_name"])

        if variant and not variant.make_vcf(ref=args["reference"], bam=bam_file, vcf_file=vcf_file, make_aux=make_aux):
            logging.error("VCF was not created.")
            return 2

        # Remove BAM file if it was generated and not kept.
        if not args["bam"] and not args["keep_temp"]:
            logging.debug("Removing BAM file: %s", bam_file)
            os.unlink(bam_file)
            os.unlink("%s.bai" % bam_file)
    else:
        vcf_file = None

    annotators_metadata = []
    if args["annotators"] and vcf_file:
        logging.info("Annotating")
        for annotator in args["annotators"]:
            # TODO: This iterates over the VCF for each annotator. Not good.
            annotator.annotate(vcf_path=vcf_file)

            meta = annotator.get_meta()

            if meta:
                annotators_metadata.append(meta)

    if args["filters"] and vcf_file:
        logging.info("Applying filters: %s", [str(f) for f in args["filters"]])
        var_set = VariantSet(vcf_file, filters=args["filters"], reference=args["reference"])

        var_set.add_metadata(mapper.get_meta())
        var_set.add_metadata(variant.get_meta())

        if args.get("version") is not None:
            var_set.add_metadata(OrderedDict({"PHEnix-Version":(args["version"],)}))

        for annotator_md in annotators_metadata:
            var_set.add_metadata(annotator_md)

        final_vcf = os.path.join(args["outdir"], "%s.filtered.vcf" % args["sample_name"])
        var_set.filter_variants()
        var_set.write_variants(final_vcf)

        if args['json'] == True:
            var_set.write_to_json(final_vcf, args["json_info"])

        if not args["vcf"] and not args["keep_temp"]:
            logging.debug("Removing VCF file: %s", vcf_file)
            os.unlink(vcf_file)

    if args["workflow"] and args["input"]:
        component_complete = os.path.join(args["outdir"], "ComponentComplete.txt")
        open(component_complete, 'a').close()

        logging.info("PIPELINE_END")

    return 0

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
