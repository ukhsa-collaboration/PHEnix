'''Implementation of the Mapper class using Bowtie2 mapper.

Created on 17 Sep 2015

@author: alex
'''
from collections import OrderedDict
import logging
import os
import subprocess
import tempfile

from phe.mapping import Mapper


class Bowtie2Mapper(Mapper):
    '''Bowtie2 mapper'''

    _default_options = "-p 1"
    """Default options for the mapper."""
    _cmd = "bowtie2"
    """Command for calling mapper."""

    name = "bowtie2"
    """Plain text name of the mapper."""

    def __init__(self, cmd_options=None):
        """Constructor for BWA mapper."""

        if not cmd_options:
            cmd_options = self._default_options

        super(Bowtie2Mapper, self).__init__(cmd_options=cmd_options)

        self.last_command = ""

    def create_aux_files(self, ref):
        if os.system("bowtie2-build %s %s" % (ref, ref)) == 0:
            return True
        else:
            return False

    def make_bam(self, *args, **kwargs):
        with tempfile.NamedTemporaryFile(suffix=".sam") as tmp:
            out_file = kwargs.get("out_file").replace(".bam", "")

            kwargs["out_file"] = tmp.name

            success = self.make_sam(*args, **kwargs)
            if not success:
                logging.warn("Could not map reads to the reference.")
                return False

            # Convert reads sam to bam filtering on MQ > 0.
            cmd = "samtools view -bhS %s | samtools sort - %s" % (tmp.name, out_file)  #  samtools view -bq 1 -
            success = os.system(cmd)
            if success != 0:
                logging.warn("Could not convert to BAM")
                logging.warn("CMD: %s", cmd)
                return False

            self.last_command += " && %s" % cmd
            cmd = "samtools index %s.bam" % out_file

            success = os.system(cmd)
            if success != 0:
                logging.warn("Could not index the BAM.")
                logging.warn("CMD: %s", cmd)
                return False
        return True


    def make_sam(self, *args, **kwargs):
        ref = kwargs.get("ref")
        r1 = kwargs.get("R1")
        r2 = kwargs.get("R2")
        out_file = kwargs.get("out_file")
        sample_name = kwargs.get("sample_name", "test_sample")

        if ref is None or r1 is None or r2 is None or out_file is None:
            logging.error("One of the required parameters is not specified.")
            return False

        d = {"cmd": self._cmd,
             "ref": os.path.abspath(ref),
             "r1": os.path.abspath(r1),
             "r2": os.path.abspath(r2),
             "out_sam": os.path.abspath(out_file),
             "sample_name": sample_name,
             "extra_options": self.cmd_options
             }

#         if self.create_aux_files(ref):
#             logging.error("Computing index has failed. Abort")
#             return False
        # TODO: should the above command have -k 1 as default option?
        cmd = "%(cmd)s --rg-id '%(sample_name)s' --rg 'SM:%(sample_name)s' %(extra_options)s -x %(ref)s -1 %(r1)s -2 %(r2)s -S %(out_sam)s" % d

        if os.system(cmd) != 0:
            logging.error("Mapping reads has failed.")

            return False

        self.last_command = cmd
        return True

    def get_version(self):

        p = subprocess.Popen(["bowtie2", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, _) = p.communicate()

        try:
            version = output.split("\n")[0].split(" ")[-1]
        except Exception:
            version = "n/a"

        return version

    def get_info(self, plain=False):

        d = {"name": "bowtie2", "version": self.get_version(), "command": self.last_command}

        if plain:
            result = "Bowtie2(%(version)s): %(command)s" % d
        else:
            result = OrderedDict(d)

        return result
