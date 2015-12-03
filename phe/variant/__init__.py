"""Classes and methods to work with variants and such."""
import abc
from collections import OrderedDict
import gzip
import logging
import pickle

from vcf import filters
import vcf
from vcf.parser import _Filter

from phe.metadata import PHEMetaData
from phe.variant_filters import make_filters, PHEFilterBase, str_to_filters


class VCFTemplate(object):
    """This is a small hack class for the Template used in generating
    VCF file."""

    def __init__(self, vcf_reader):
        self.infos = vcf_reader.infos
        self.formats = vcf_reader.formats
        self.filters = vcf_reader.filters
        self.alts = vcf_reader.alts
        self.contigs = vcf_reader.contigs
        self.metadata = vcf_reader.metadata
        self._column_headers = vcf_reader._column_headers
        self.samples = vcf_reader.samples

class VariantSet(object):
    """A convenient representation of set of variants.
    TODO: Implement iterator and generator for the variant set.
    """

    _reader = None

    def __init__(self, vcf_in, filters=None):
        """Constructor of variant set.
        
        Parameters:
        -----------
        vcf_in: str
            Path to the VCF file for loading information.
        filters: str or dict, optional
            Dictionary or string of the filter:threshold key value pairs.
        """
        self.vcf_in = vcf_in
        self._reader = vcf.Reader(filename=vcf_in)
        self.out_template = VCFTemplate(self._reader)

        self.filters = []
        if filters is not None:
            if isinstance(filters, str):
                self.filters = str_to_filters(filters)
            elif isinstance(filters, dict):
                self.filters = make_filters(config=filters)
            elif isinstance(filters, list):
                self.filters = filters
            else:
                logging.warn("Could not create filters from %s", filters)
        else:
            reader = vcf.Reader(filename=self.vcf_in)
            filters = {}
            for filter_id in reader.filters:
                filters.update(PHEFilterBase.decode(filter_id))

            if filters:
                self.filters = make_filters(config=filters)

        self.variants = []

    def filter_variants(self, keep_only_snps=True):
        """Create a variant """

        if self._reader is None:
            # Create a reader class from input VCF.
            self._reader = vcf.Reader(filename=self.vcf_in)

        # get list of existing filters.
        existing_filters = {}
        removed_filters = []

        for filter_id in self._reader.filters:
            conf = PHEFilterBase.decode(filter_id)
            tuple(conf.keys())
            existing_filters.update({tuple(conf.keys()):filter_id})

        # Add each filter we are going to use to the record.
        # This is needed for writing out proper #FILTER header in VCF.
        for record_filter in self.filters:
            # We know that each filter has short description method.
            short_doc = record_filter.short_desc()
            short_doc = short_doc.split('\n')[0].lstrip()

            filter_name = PHEFilterBase.decode(record_filter.filter_name())

            # Check if the sample has been filtered for this type of filter
            #    in the past. If so remove is, because it is going to be refiltered.
            if tuple(filter_name) in existing_filters:
                logging.info("Removing existing filter: %s", existing_filters[tuple(filter_name)])
                removed_filters.append(existing_filters[tuple(filter_name)])
                del self._reader.filters[existing_filters[tuple(filter_name)]]

            self._reader.filters[record_filter.filter_name()] = _Filter(record_filter.filter_name(), short_doc)

        # For each record (POSITION) apply set of filters.
        for record in self._reader:

            # If this record failed filters and we removed some,
            #    check is they need to be removed from record.
            if isinstance(record.FILTER, list) and len(record.FILTER) > 0:
                for filter_id in removed_filters:
                    if filter_id in record.FILTER:
                        record.FILTER.remove(filter_id)

            for record_filter in self.filters:

                # Call to __call__ method in each filter.
                result = record_filter(record)

                # Record is KEPT if filter returns None
                if result == None:
                    continue

                # If we got this far, then record is filtered OUT.
                record.add_filter(record_filter.filter_name())

            # After applying all filters, check if FILTER is None.
            # If it is, then record PASSED all filters.
            if record.FILTER is None or record.FILTER == []:
                record.FILTER = 'PASS'
                # FIXME: Does this work for indels?
                if keep_only_snps and record.is_snp:
                    self.variants.append(record)
            else:
                self.variants.append(record)

        self.update_filters(self._reader.filters)

        return [ variant for variant in self.variants if variant.FILTER == "PASS"]

    def add_metadata(self, info):
        """Add metadata to the variant set.
        
        Parameters:
        -----------
        info: dict
            Dictionary of key value pairs to be inserted into metadata.
        """
        for info_key, metadata in info.items():
            self.out_template.metadata[info_key] = metadata

    def write_variants(self, vcf_out, only_snps=False, only_good=False):
        """Write variants to a VCF file.
        
        Parameters:
        -----------
        vcf_out: str
            Path to the file where VCF data is written.
        only_snps: bool, optional
            True is *only* SNP are to be written to the output (default: False).
        only_good: bool, optional
            True if only those records that PASS all filters should be written
            (default: False).
        
        Returns:
        int:
            Number of records written.
        """
        written_variants = 0

        # Check if the output file ends with .gz, then compress data.
        open_func = gzip.open if vcf_out.endswith(".gz") else open

        with open_func(vcf_out, "w") as out_vcf:
            writer = vcf.Writer(out_vcf, self.out_template)

            # Output internal variants (if exist) otherwise, output data from reader.
            variants = self.variants if self.variants else self._reader

            for record in variants:

                if only_snps and not record.is_snp:
                    continue

                if only_good and record.FILTER:
                    continue

                writer.write_record(record)
                written_variants += 1

        return written_variants

    def _write_bad_variants(self, vcf_out):
        """**PRIVATE:** Write only those records that **haven't** passed."""
        written_variants = 0
        # Check if the output file ends with .gz, then compress data.
        open_func = gzip.open if vcf_out.endswith(".gz") else open
        with open_func(vcf_out, "w") as out_vcf:
            writer = vcf.Writer(out_vcf, self.out_template)
            for record in self.variants:
                if record.FILTER != "PASS" and record.FILTER is not None:
                    writer.write_record(record)
                    written_variants += 1
        return written_variants

    def serialise(self, out_file):
        """Save the data in this class to a file for future use/reload.
        
        Parameters:
        -----------
        out_file: str
            path to file where the data should be written to.
            
        Returns:
        --------
        int:
            Number of variants written.
        """
        written_variants = 0
        open_func = gzip.open if out_file.endswith(".gz") else open
        with open_func(out_file, "w") as out_vcf:
            writer = vcf.Writer(out_vcf, self.out_template)
            for record in self.variants:
                writer.write_record(record)
                written_variants += 1

        return written_variants

    def update_filters(self, new_filters):
        """Update internal filters in the output template."""
        for new_filter, filter_data in new_filters.items():
            self.out_template.filters[new_filter] = filter_data


class VariantCaller(PHEMetaData):
    """Abstract class used for access to the implemented variant callers."""

    __metaclass__ = abc.ABCMeta

    def __init__(self, cmd_options=None):
        """Constructor for variant caller.
        
        Parameters:
        -----------
        cmd_options: str, optional
            Command options to pass to the variant command.
        """
        self.cmd_options = cmd_options

        super(VariantCaller, self).__init__()

    @abc.abstractmethod
    def make_vcf(self, *args, **kwargs):
        """Create a VCF from **BAM** file.
        
        Parameters:
        -----------
        ref: str
            Path to the reference file.
        bam: str
            Path to the indexed **BAM** file for calling variants.
        vcf_file: str
            path to the VCF file where data will be written to.
            
        Returns:
        --------
        bool:
            True if variant calling was successful, False otherwise.
        """
        raise NotImplementedError("make_vcf is not implemented yet.")

    @abc.abstractmethod
    def create_aux_files(self, ref):
        """Create needed (if any) auxiliary files.
        These files are required for proper functioning of the variant caller.
        """
        raise NotImplementedError("create_aux_files is not implemeted.")

    @abc.abstractmethod
    def get_info(self, plain=False):
        """Get information about this variant caller."""
        raise NotImplementedError("Get info has not been implemented yet."
                                  )
    def get_meta(self):
        """Get the metadata about this variant caller."""
        od = self.get_info()
        od["ID"] = "VariantCaller"
        return OrderedDict({"PHEVariantMetaData": [od]})

    @abc.abstractmethod
    def get_version(self):
        """Get the version of the underlying command used."""
        raise NotImplementedError("Get version has not been implemented yet.")
