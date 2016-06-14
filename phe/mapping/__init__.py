"""Mapping related classes and functions.

.. _implementing-mapper:

Implementing Mapper
-------------------

This project does not provide novel way of mapping or calling variants. Instead
it interfaces with available mappers through command line interface. Hopefully,
we have made it straight forward to add your favourite mapper for others to use.

To do this, you will need to implement :py:class:`Mapper` class. There are a
number of methods that need to be implemented for you mapper to work:

- :py:attr:`Mapper.name` - Attribute that specifies the name of the mapper (used to dynamically load it).
- :py:meth:`Mapper.make_sam` - Create a sam file in *out_file* from *R1/Forward*, *R2/Reverse* and *ref*. See :py:meth:`Mapper.make_sam` for all input options.
- :py:meth:`Mapper.create_aux_files` - Create auxilliary files needed for the mapper. This is also called by :ref:`prepare-script`.
- :py:meth:`Mapper.get_version` - get the version of the mapper.
- :py:meth:`Mapper.get_info` - Get the meta data information for the mapper (included in the VCF header).
- :py:meth:`Mapper.make_bam` - If you want, you can also implement this method, but my default samtools is used to convert from SAM -> BAM.

Bear in mind that **cmd_options** will be passed on, if they specified in the command line or config under **mapper-options**.

Once you have implemented this interface, save the file in the *mapping*
directory and it should be automatically picked up. To verify run:

.. code-block:: bash

   run_snp_pipeline.py --help
   
You should see your mapper in the list of available mappers.
If the mapper works, it will also be included n automatically
in the documentations.


:Date: 22 Sep, 2015
:Author: Alex Jironkin
"""

import abc
from collections import OrderedDict
import logging
import os
import shlex
from subprocess import Popen
import subprocess
import tempfile

from phe.metadata import PHEMetaData
from phe.utils import calculate_memory_for_sort


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
        
        Parameters
        ----------
        cmd_options: str, optional
            Command options for the mapper.
        """
        super(Mapper, self).__init__()

        self.cmd_options = cmd_options

        self.validate()

    @abc.abstractmethod
    def create_aux_files(self, ref):
        """Create required auxiliary files for *ref*.
        
        Parameters
        ----------
        ref: str
            Path to the reference file.
        
        Returns
        -------
        bool:
            True if auxiliary files were created, False otherwise.
        """
        raise NotImplementedError("creating auxiliary has not been implemented yet.")

    @abc.abstractmethod
    def make_sam(self, *args, **kwargs):
        """Make SAM from reference, and fastq files.
        
        Parameters
        ----------
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

        Returns
        -------
        bool:
            True iff mapping returns 0, False otherwise.
        """
        raise NotImplementedError("make_sam is not implemented yet.")

    def make_bam(self, *args, **kwargs):
        """Make **indexed** BAM from reference, and fastq files.

        Parameters
        ----------
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

        Returns
        -------
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
            sort_memory = calculate_memory_for_sort()


            with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:

                if samtools_version[0] >= 1 and samtools_version[1] >= 3:
                    out_file += ".bam"
                    view_cmd = "samtools view -bhS %s" % tmp.name  #  samtools view -bq 1 -
                    if sort_memory is None:
                        sort_cmd = "samtools sort %s -o %s" % (tmp_bam.name, out_file)
                    else:
                        sort_cmd = "samtools sort -m %s %s -o %s" % (sort_memory, tmp_bam.name, out_file)
                else:
                    view_cmd = "samtools view -bhS %s" % tmp.name  #  samtools view -bq 1 -
                    if sort_memory is None:
                        sort_cmd = "samtools sort %s %s" % (tmp_bam.name, out_file)
                    else:
                        sort_cmd = "samtools sort -m %s %s %s" % (sort_memory, tmp_bam.name, out_file)


                logging.debug("SAMTOOLS VERSION: %s, CMD: %s && %s", samtools_version, view_cmd, sort_cmd)

                p = Popen(shlex.split(view_cmd), stderr=subprocess.PIPE, stdout=tmp_bam)
                (_, view_stderr) = p.communicate()

                if p.returncode != 0:
                    logging.error("Could not convert to BAM: %s", view_cmd)
                    return False

                p = Popen(shlex.split(sort_cmd), stderr=subprocess.PIPE, stdout=subprocess.PIPE)


                (stdout, sort_stderr) = p.communicate()

                if p.returncode != 0:

                    logging.error("Could not convert to BAM: %s", sort_cmd)
                    logging.error("%s", sort_stderr)
                    return False

            self.last_command += " && %s && %s" % (view_cmd, sort_cmd)
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


    def validate(self):
        """Validate itself, by checking appropriate commands can run."""

        if self.get_version() == "n/a":
            raise Exception("%s is not available in your PATH." % self.name)
