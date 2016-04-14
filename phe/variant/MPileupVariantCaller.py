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

        if p.returncode != 0:
            version = "n/a"
        else:
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

        try:
            version = [ int(v) for v in self.get_version().split(".")]
            if len(version) == 2:
                version.append(0)
        except ValueError:
            # Older versions of samtools don't have --version command
            version = [0, 0, 0]

        with tempfile.NamedTemporaryFile(suffix=".pileup") as tmp:
            opts["pileup_file"] = tmp.name

            if version[0] >= 1 and version[1] >= 3:
                pileup_cmd = "samtools mpileup -t DP,AD,SP,ADF,ADR,INFO/AD,INFO/ADF,INFO/ADR -Auf %(ref)s -o %(pileup_file)s %(bam)s" % opts
            else:
                pileup_cmd = "samtools mpileup -t DP,DV,DP4,DPR,SP -Auf %(ref)s -o %(pileup_file)s %(bam)s" % opts

            bcf_cmd = "bcftools call %(extra_cmd_options)s -o %(all_variants_file)s %(pileup_file)s" % opts

            logging.debug("CMD: %s", pileup_cmd)
            p = Popen(shlex.split(pileup_cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # FIXME: Strictly speaking those should be p.stderr.readline(), but that reads 1 char at a time.
            stderr = []
            for line in p.stderr:
                line = line.strip()
                logging.debug(line)
                stderr.append(line)

            p.wait()

            if p.returncode != 0:
                logging.error("Pileup creation failed.")
                logging.error("\n".join(stderr))
                return False

            p = Popen(shlex.split(bcf_cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            stderr = []
            for line in p.stderr:
                line = line.strip()
                logging.debug(line)
                stderr.append(line)

            p.wait()

            if p.returncode != 0:
                logging.warn("Pileup VCF creation was not successful.")
                logging.error("\n".join(stderr))
                return False

            self.last_command = "%s && %s" % (pileup_cmd, bcf_cmd)

        return True

    def create_aux_files(self, ref):
        """Index reference with faidx from samtools.

        Parameters
        ----------
        ref: str
            Path to the reference file.

        Returns
        -------
        bool:
            True if auxiliary files were created, False otherwise.
        """

        p = Popen(["samtools", "faidx", ref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if p.returncode != 0:
            logging.warn("Fasta index could not be created.")
            return False

        return True
