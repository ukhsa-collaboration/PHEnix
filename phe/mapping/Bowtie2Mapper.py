'''Implementation of the Mapper class using Bowtie2 mapper.

:Date: 17 Sep, 2015
:Author: Alex Jironkin
'''
from collections import OrderedDict
import logging
import os
import shlex
from subprocess import Popen
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

        self.last_command = None

    def create_aux_files(self, ref):
        p = Popen(["bowtie2-build", ref, ref], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()

        if p.returncode == 0:
            return True
        else:
            return False

    def make_sam(self, *args, **kwargs):
        ref = kwargs.get("ref")
        r1 = kwargs.get("R1")
        r2 = kwargs.get("R2")
        out_file = kwargs.get("out_file")
        sample_name = kwargs.get("sample_name", "test_sample")

        make_aux = kwargs.get("make_aux", False)

        if ref is None or r1 is None or r2 is None or out_file is None:
            logging.error("One of the required parameters is not specified.")
            return False

        d = {"cmd": self._cmd,
             "ref": os.path.abspath(ref),
             "r1": os.path.abspath(r1),
             "r2": os.path.abspath(r2),
             "out_sam": out_file,
             "sample_name": sample_name,
             "extra_options": self.cmd_options
             }

        if make_aux:
            if not self.create_aux_files(ref):
                logging.error("Computing index has failed. Abort")
                return False

        # TODO: should the above command have -k 1 as default option?
        cmd = r"%(cmd)s --rg-id '%(sample_name)s' --rg 'SM:%(sample_name)s' %(extra_options)s -x %(ref)s -1 %(r1)s -2 %(r2)s" % d

        logging.debug("CMD: %s", cmd)

        p = Popen(shlex.split(cmd), stdout=d["out_sam"], stderr=subprocess.PIPE)

        stderr = []
        for line in p.stderr:
            line = line.strip()
            logging.debug(line)
            stderr.append(line)

        p.wait()

        if p.returncode != 0:
            logging.error("Mapping reads has failed.")
            logging.error("\n".join(stderr))

            return False

        self.last_command = cmd
        return True

    def get_version(self):
        try:
            p = subprocess.Popen(["bowtie2", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            (output, _) = p.communicate()
        except OSError as e:
            logging.error(str(e))
            return "n/a"

        if p.returncode != 0:
            version = "n/a"
        else:
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
