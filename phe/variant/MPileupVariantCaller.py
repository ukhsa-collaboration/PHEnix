'''
Created on 22 Sep 2015

@author: alex
'''
from collections import OrderedDict
import logging
import os
import subprocess
import tempfile

from phe.variant import VariantCaller


class MPileupVariantCaller(VariantCaller):
    """Implemetation of the Broad institute's variant caller."""

    name = "mpileup"
    """Plain text name of the variant caller."""

    _default_options = "-m -f GQ"
    """Default options for the variant caller."""

    def __init__(self, cmd_options=None):
        """Constructor"""
        if cmd_options is None:
            cmd_options = self._default_options

        super(MPileupVariantCaller, self).__init__(cmd_options=cmd_options)

        self.last_command = None

    def get_info(self, plain=False):
        d = {"name": self.name, "version": self.get_version(), "command": self.last_command}

        if plain:
            result = "mpileup(%(version)s): %(command)s" % d
        else:
            result = OrderedDict(d)

        return result

    def get_version(self):

        p = subprocess.Popen(["samtools", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, _) = p.communicate()

        # first line is the version of the samtools
        version = output.split("\n")[0].split(" ")[1]

        return version

    def make_vcf(self, *args, **kwargs):
        ref = kwargs.get("ref")
        bam = kwargs.get("bam")

        make_aux = kwargs.get("make_aux", False)

        if kwargs.get("vcf_file") is None:
            kwargs["vcf_file"] = "variants.vcf"

        opts = {"ref": os.path.abspath(ref),
                "bam": os.path.abspath(bam),
                "all_variants_file": os.path.abspath(kwargs.get("vcf_file")),
                "extra_cmd_options": self.cmd_options}

        if make_aux:
            if not self.create_aux_files(ref):
                logging.warn("Auxiliary files were not created.")
                return False

        with tempfile.NamedTemporaryFile(suffix=".pileup") as tmp:
            opts["pileup_file"] = tmp.name
            cmd = "samtools mpileup -t DP,DV,DP4,DPR,SP -Auf %(ref)s %(bam)s | bcftools call %(extra_cmd_options)s > %(all_variants_file)s" % opts
            print cmd
            self.last_command = cmd
            if os.system(cmd) != 0:
                logging.warn("Pileup creation was not successful.")
                return False

        return True

    def create_aux_files(self, ref):
        """Index reference with faidx from samtools.

        Parameters:
        -----------
        ref: str
            Path to the reference file.

        Returns:
        --------
        bool:
            True if auxiliary files were created, False otherwise.
        """

        success = os.system("samtools faidx %s" % ref)

        if success != 0:
            logging.warn("Fasta index could not be created.")
            return False

        return True
