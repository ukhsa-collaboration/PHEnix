"""Mapping related classes and functions."""

import abc
from collections import OrderedDict

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

        Returns:
        --------
        bool:
            True iff mapping returns 0, False otherwise.
        """
        raise NotImplementedError("make_sam is not implemented yet.")

    @abc.abstractmethod
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
        raise NotImplementedError("make_bam is not implemented yet.")

    @abc.abstractmethod
    def get_info(self, plain=False):
        """Get information about the mapper."""
        raise NotImplementedError("get_info is not implemented yet.")

    def get_meta(self):
        od = self.get_info()
        od["ID"] = "Mapper"
        return OrderedDict({"PHEMapperMetaData": [od]})

    @abc.abstractmethod
    def get_version(self):
        """Get the version of the underlying command used."""
        raise NotImplementedError("Get version has not been implemented yet.")
