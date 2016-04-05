"""Mapping related classes and functions."""

import abc
from collections import OrderedDict
import logging
import os
import shlex
from subprocess import Popen
import subprocess
import tempfile

from phe.metadata import PHEMetaData


class Mapper(PHEMetaData):
    """Abstract Mapper class that provides generic interface to the 
    mapping implementation for a particular mapper.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        """Return name of the mapper (implemented by individual classes."""
        return self.name

    def __init__(self, cmd_options=None):
        """Abstract constructor for the Mapper class.
        
        Parameters:
        -----------
        cmd_options: str, optional
            Command options for the mapper.
        """
        super(Mapper, self).__init__()

        self.cmd_options = cmd_options

    @abc.abstractmethod
    def create_aux_files(self, ref):
        """Create required auxiliary files for *ref*.

        Parameters:
        -----------
        ref: str
            Path to the reference file.
        
        Returns:
        --------
        bool:
            True if auxiliary files were created, False otherwise.
        """
        raise NotImplementedError("creating auxiliary has not been implemented yet.")

    @abc.abstractmethod
    def make_sam(self, *args, **kwargs):
        """Make SAM from reference, and fastq files.

        Parameters:
        -----------
        ref: str
            Path to the reference file used for mapping.
        R1: str
            Path to the R1/Forward reads file in FastQ format.
        R2: str:
            Path to the R2/Reverse reads file in FastQ format.
        out_file: str
            Name of the output file where SAM data is written.
        sample_name: str
            Name of the sample to be included in the read group (RG) field.
        make_aux: bool, optional
            True/False to make auxilliary files (default: False).

        Returns:
        --------
        bool:
            True iff mapping returns 0, False otherwise.
        """
        raise NotImplementedError("make_sam is not implemented yet.")

    def make_bam(self, *args, **kwargs):
        """Make **indexed** BAM from reference, and fastq files.

        Parameters:
        -----------
        ref: str
            Path to the reference file used for mapping.
        R1: str
            Path to the R1/Forward reads file in FastQ format.
        R2: str:
            Path to the R2/Reverse reads file in FastQ format.
        out_file: str
            Name of the output file where BAM data is written.
        sample_name: str
            Name of the sample to be included in the read group (RG) field.

        Returns:
        --------
        bool:
            True iff mapping, converting to BAM and indexing returns 0,
            False otherwise.
        """
        with tempfile.NamedTemporaryFile(suffix=".sam") as tmp:
            out_file = kwargs.get("out_file").replace(".bam", "")

            kwargs["out_file"] = tmp

            success = self.make_sam(*args, **kwargs)
            if not success:
                logging.warn("Could not map reads to the reference.")
                return False

            # Convert reads sam to bam filtering on MQ > 0.
            samtools_version = self.get_samtools_version()

            if samtools_version[0] >= 1 and samtools_version[1] >= 3:
                out_file += ".bam"
                view_cmd = "samtools view -bhS %s" % tmp.name  #  samtools view -bq 1 -
                sort_cmd = "samtools sort - -o %s" % out_file
            else:
                view_cmd = "samtools view -bhS %s" % tmp.name  #  samtools view -bq 1 -
                sort_cmd = "samtools sort - %s" % out_file

            logging.debug("SAMTOOLS VERSION: %s, CMD: %s | %s", samtools_version, view_cmd, sort_cmd)

            # In order to change to Popen need to split the pipe into 2 processes.
            p = {"view":None, "sort":None}

            p["view"] = Popen(shlex.split(view_cmd), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            p["sort"] = Popen(shlex.split(sort_cmd), stdin=p["view"].stdout, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

            (_, view_stderr) = p["view"].communicate()
            (stdout, sort_stderr) = p["sort"].communicate()

            if p["view"].returncode != 0 or p["sort"].returncode != 0:

                logging.error("Could not convert to BAM: %s | %s", shlex.split(view_cmd), shlex.split(sort_cmd))
                logging.error("%s", view_stderr)
                logging.error("%s", sort_stderr)
                return False

            self.last_command += " && %s | %s" % (view_cmd, sort_cmd)
            if not out_file.endswith(".bam"):
                out_file += ".bam"

            cmd = ["samtools", "index", out_file]
            p = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdout, stderr) = p.communicate()

            if p.returncode != 0:
                logging.warn("Could not index the BAM.")
                logging.warn("CMD: %s", " ".join(cmd))
                logging.warn(stderr)
                return False
        return True

    @abc.abstractmethod
    def get_info(self, plain=False):
        """Get information about the mapper."""
        raise NotImplementedError("get_info is not implemented yet.")

    def get_meta(self):
        od = self.get_info()
        od["ID"] = "Mapper"
        return OrderedDict({"PHEMapperMetaData": [od]})

    def get_samtools_version(self):
        """Get version of samtools used. Reterned as a triple (major, minor, patch)"""

        p = subprocess.Popen(["samtools", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, _) = p.communicate()

        # first line is the version of the samtools

        try:
            version = [ int(v) for v in output.split("\n")[0].split(" ")[1].split(".")]
            if len(version) == 2:
                version.append(0)
        except ValueError:
            # Older versions of samtools don't have --version command
            version = [0, 0, 0]

        return tuple(version)

    @abc.abstractmethod
    def get_version(self):
        """Get the version of the underlying command used."""
        raise NotImplementedError("Get version has not been implemented yet.")
