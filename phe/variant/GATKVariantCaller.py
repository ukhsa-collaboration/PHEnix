'''
:Date: 22 Sep, 2015
:Author: Alex Jironkin
'''
from collections import OrderedDict
import logging
import os
import shlex
from subprocess import Popen
import subprocess

from phe.variant import VariantCaller


class GATKVariantCaller(VariantCaller):
    """Implemetation of the Broad institute's variant caller."""

    name = "gatk"
    """Plain text name of the variant caller."""

    _default_options = "--sample_ploidy 2 --genotype_likelihoods_model SNP -rf BadCigar -out_mode EMIT_ALL_SITES -nt 1"
    """Default options for the variant caller."""

    def __init__(self, cmd_options=None):
        """Constructor"""
        if cmd_options is None:
            cmd_options = self._default_options

        super(GATKVariantCaller, self).__init__(cmd_options=cmd_options)

        self.last_command = None

        self.validate()

    def get_info(self, plain=False):
        d = {"name": "gatk", "version": self.get_version(), "command": self.last_command}

        if plain:
            result = "GATK(%(version)s): %(command)s" % d
        else:
            result = OrderedDict(d)

        return result

    def get_version(self):

        p = subprocess.Popen(["java", "-jar", os.environ["GATK_JAR"], "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, _) = p.communicate()

        if p.returncode != 0:
            version = "n/a"
        else:
            # last character is EOL.
            version = output.split("\n")[-2]

        return version

    def make_vcf(self, *args, **kwargs):
        ref = kwargs.get("ref")
        bam = kwargs.get("bam")

        make_aux = kwargs.get("make_aux", False)

        if kwargs.get("vcf_file") is None:
            kwargs["vcf_file"] = "variants.vcf"

        opts = {"ref": os.path.abspath(ref),
                "bam": os.path.abspath(bam),
                "gatk_jar": os.environ["GATK_JAR"],
                "all_variants_file": os.path.abspath(kwargs.get("vcf_file")),
                "extra_cmd_options": self.cmd_options}
        if make_aux:
            if not self.create_aux_files(ref):
                logging.warn("Auxiliary files were not created.")
                return False

        # Call variants
        # FIXME: Sample ploidy = 2?
        cmd = "java -XX:+UseSerialGC -jar %(gatk_jar)s -T UnifiedGenotyper -R %(ref)s -I %(bam)s -o %(all_variants_file)s %(extra_cmd_options)s" % opts
        logging.debug("CMD: %s", cmd)

        p = Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        (stdout, stderr) = p.communicate()

        for line in stdout.split("\n"):
            logging.debug(line)
        for line in stderr.split("\n"):
            logging.debug(line)

        if p.returncode != 0:
            logging.warn("Calling variants returned non-zero exit status.")
            logging.warn(stderr)
            return False

        self.last_command = cmd

        return True

    def create_aux_files(self, ref):
        """Create auxiliary files needed for this variant.

        Tools needed: samtools and picard tools. Picard tools is a Java
        library that can be defined using environment variable: PICARD_JAR
        specifying path to picard.jar or PICARD_TOOLS_PATH specifying path
        to the directory where separate jars are (older version before jars
        were merged into a single picard.jar).
        
        Parameters
        ----------
        ref: str
            Path to the reference file.
        
        Returns
        -------
        bool:
            True if auxiliary files were created, False otherwise.
        """

        ref_name, _ = os.path.splitext(ref)

        p = Popen(["samtools", "faidx", ref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()

        if p.returncode != 0:
            logging.warn("Fasta index could not be created.")
            logging.warn(stderr)
            return False

        d = {"ref": ref, "ref_name": ref_name}

        if os.environ.get("PICARD_TOOLS_PATH"):
            d["picard_tools_path"] = os.path.join(os.environ["PICARD_TOOLS_PATH"], "CreateSequenceDictionary.jar")
        elif os.environ.get("PICARD_JAR"):
            # This is used in newer version of PICARD tool where multiple
            #    jars were merged into a single jar file.
            d["picard_tools_path"] = "%s %s" % (os.environ["PICARD_JAR"], "CreateSequenceDictionary")
        else:
            logging.error("Picard tools are not present in the path. Please set PICARD_TOOLS_PATH to CreateSequenceDictionary.jar for an older many commands version or PICARD_JAR for single jar file.")
            return False

        if not os.path.exists("%s.dict" % ref_name):
            cmd = "java -jar %(picard_tools_path)s R=%(ref)s O=%(ref_name)s.dict" % d
            p = Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdout, stderr) = p.communicate()
            if p.returncode != 0:
                logging.warn("Dictionary for the %s reference could not be created", ref)
                logging.warn(stderr)
                return False
        else:
            logging.debug("PICARD AUX EXISTS.")

        return True

    def validate(self):
        if "GATK_JAR" not in os.environ:
            raise Exception("GATK_JAR is not found in you environment. Please set environemnt variable(s) and run the tool again.")
