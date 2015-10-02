'''
Created on 17 Sep 2015

@author: alex
'''
import logging
import os
import sys
import tempfile

from phe.mapping import Mapper


class BWAMapper(Mapper):
    '''
    classdocs
    '''

    _default_options = "-t 1"
    _cmd = "bwa mem"

    name = "bwa"

    def __init__(self, cmd_options=None):
        '''
        Constructor
        '''

        if not cmd_options:
            cmd_options = self._default_options

        super(BWAMapper, self).__init__(cmd_options=cmd_options)

    def make_bam(self, *args, **kwargs):
        with tempfile.NamedTemporaryFile(suffix=".sam") as tmp:
            out_file = kwargs.get("out_file").replace(".bam", "")

            kwargs["out_file"] = tmp.name

            success = self.make_sam(*args, **kwargs)
            if not success:
                logging.warn("Could not map reads to the reference.")
                return False

            cmd = "samtools view -bhS %s | samtools view -bq 2 - | samtools sort - %s" % (tmp.name, out_file)
            success = os.system(cmd)
            if success != 0:
                logging.warn("Could not convert to BAM")
                logging.warn("CMD: %s", cmd)
                return False

            cmd = "samtools index %s.bam" % out_file
            success = os.system(cmd)
            if success != 0:
                logging.warn("Could not index the BAM.")
                logging.warn("CMD: %s", cmd)
                return False
        return True


    def make_sam(self, *args, **kwargs):

        ref = kwargs.get("ref")
        R1 = kwargs.get("R1")
        R2 = kwargs.get("R2")
        out_file = kwargs.get("out_file")
        sample_name = kwargs.get("sample_name", "test_sample")

        if ref is None or R1 is None or R2 is None or out_file is None:
            logging.error("One of the required parameters is not specified.")
            return False

        d = {"cmd": self._cmd,
             "ref": ref,
             "r1": R1,
             "r2": R2,
             "out_sam": out_file,
             "sample_name": sample_name,
             "extra_options": self.cmd_options
             }

        if os.system("bwa index %(ref)s" % d) != 0:
            logging.error("Computing index has failed. Abort")
            return False

        cmd = "%(cmd)s -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s' %(extra_options)s %(ref)s %(r1)s %(r2)s > %(out_sam)s" % d

        if os.system(cmd) != 0:
            logging.error("Mapping reads has failed.")
            return False

        return True

    def get_info(self):
        return None
