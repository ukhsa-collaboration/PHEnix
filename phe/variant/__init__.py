"""Classes and methods to work with _variants and such.

.. _implementing-variant:

Implementing Variant Caller
---------------------------


Similarly to :ref:`implementing-mapper` and :ref:`implementing-filter`
adding your variant caller is easy and sraight forward. Implement
:py:class:`VariantCaller` class with the following attributes/methods:

- :py:attr:`VariantCaller.name` - Name of your variant caller, used to dynamically load it.
- :py:meth:`VariantCaller.make_vcf` - Create a VCF from input BAM. See :py:meth:`VariantCaller.make_vcf` for input list.
- :py:meth:`VariantCaller.get_version` - Return the version of this variant caller.
- :py:meth:`VariantCaller.get_info` - Return information about this caller and included in the VCF header.
- :py:meth:`VariantCaller.create_aux_files` - Create auxilliary files needed by the caller. Also used in  :ref:`prepare-script`.

Save your new variant caller in the *variant* directiry and it will be
automatically picked up by the system. To verify run:

.. code-block:: bash

   run_snp_pipeline.py --help

It should appear in the list of available callers.

:Date: 22 Sep, 2015
:Modified: |today|
:Author: Alex Jironkin
"""
import abc
from collections import OrderedDict, namedtuple
import gzip
import logging
import operator
import pickle
import json

from Bio import SeqIO
from vcf import filters
import vcf
from vcf.model import make_calldata_tuple
from vcf.parser import _Filter

from phe.metadata import PHEMetaData
from phe.variant_filters import make_filters, PHEFilterBase, str_to_filters
from phe.utils import is_uncallable


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
    """A convenient representation of set of _variants.
    TODO: Implement iterator and generator for the variant set.
    """

    _reader = None

    def __init__(self, vcf_in, filters=None, reference=None):
        """Constructor of variant set.

        Parameters
        ----------
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

        self._variants = []

        self._read_reference(reference)

    def __iter__(self):
        '''Iterator over **all** variants'''
        return self.variants(only_good=False)

    def variants(self, only_good=False):
        '''Generator over the variant set.

        Parameters
        ----------
        only_good: bool, optional
            Iff True good and bad variants are returned,
            otherwise only good are returned (default: False).
        '''
        for var in self._variants:
            if not only_good:
                yield var
            elif not var.FILTER :
                yield var

    def filter_variants(self, keep_only_snps=False, only_good=False):
        """Filter the VCF records.

        Parameters
        ----------
        keep_only_snps: bool, optional
            Retain only SNP variants (default: False).
        only_good: bool, optional
            True/False if only SNPs that PASS should output.

        Returns
        -------
        list of records is returned.
         """

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

        # Update the filters for output.
        self._update_filters(self._reader.filters)

        _pos = 1
        _chrom = None
        # For each record (POSITION) apply set of filters.
        for record in self._reader:

            if _chrom != record.CHROM:
                _pos, _chrom = 1, record.CHROM

            # Fill in any missing consecutive data with GT=./. records.
            while _pos <= record.POS:
                if _pos == record.POS:
                    _record = record
                else:
                    # This is a padding "N" record when records do not follow each other,
                    #    and there is a gap. e,g, 1,2,3,5,6 -> in 4 "N" will be inserted.

                    _ref = self._get_reference_base(record.CHROM, _pos)

                    _record = vcf.model._Record(record.CHROM, _pos, ".", _ref, [None], 0, [], {}, 'GT', None)
                    _calls = []
                    sorted_samples = sorted(record._sample_indexes.items(), key=operator.itemgetter(1))
                    for sample, i in sorted_samples:

                        _data = make_calldata_tuple(["GT"])
                        _data._types = ["String"]
                        _data._nums = [1]
                        d = ["./."]
                        _calls.append(vcf.model._Call(_record, sample=sample, data=_data(*d)))

                    _record.samples = _calls
                    _record._sample_indexes = dict(sorted_samples)


                self._filter_record(_record, removed_filters)

                # After applying all filters, check if FILTER is None.
                # If it is, then record PASSED all filters.
                if _record.FILTER is None or _record.FILTER == []:
                    if not record.is_monomorphic:
                        _record.FILTER = []

                        if not keep_only_snps or (_record.is_snp and keep_only_snps):

                            self._variants.append(_record)

                elif not only_good:
                    self._variants.append(_record)

                _pos += 1
                if _chrom is None:
                    _chrom = record.CHROM
        return [ variant for variant in self._variants if not variant.FILTER]


# --------------------------------------------------------------------------------------------------

    def _filter_record(self, record, removed_filters=list()):
        '''**PRIVATE** Filter record.

        Parameters
        ----------
        record: :py:class:vcf.Record
            Record to filter.
        '''

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
            if result is None:
                continue

            # If we got this far, then record is filtered OUT.
            record.add_filter(record_filter.filter_name())

    def _get_reference_base(self, chrom, pos):
        """Get the reference base at chromosome *chrom* and position
        *pos*. If it is not found, **N** is returned.

        Parameters
        ----------
        reference: dict
            Dictionary of chromosome to sequence map.
        chrom: str
            Chromosome to access.
        pos: int
            Position to retrieve.

        Returns
        -------
        _ref: str
            Reference base to use.
        """

        if self._reference is None:
            _ref = "N"
        else:
            try:
                _chrom = self._reference.get(chrom, [])

                # VCFs are 1-base indexed, strings are 0-base indexed.
                _ref = _chrom[pos - 1]
            except ValueError:
                logging.error("Could not retrieve reference base for %s @ %s", chrom, pos - 1)
                _ref = "N"

        return _ref

    def _read_reference(self, reference):
        """Read in the reference from the file into a dictionary.

        Parameters
        ----------
        reference: str
            Path to the reference.

        Returns
        -------
        reference: dict
            Dictionary with chromosome - sequence mapping.
        """

        if reference is None:
            self._reference = None
        else:
            self._reference = {}

            with open(reference) as reference_fp:
                for record in SeqIO.parse(reference_fp, "fasta"):
                    self._reference[record.id] = list(record.seq)

    def add_metadata(self, info):
        """Add metadata to the variant set.

        Parameters
        ----------
        info: dict
            Dictionary of key value pairs to be inserted into metadata.
        """
        for info_key, metadata in info.items():
            self.out_template.metadata[info_key] = metadata

    def write_variants(self, vcf_out, only_snps=False, only_good=False):
        """Write _variants to a VCF file.

        Parameters
        ----------
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

            # Output internal _variants (if exist) otherwise, output data from reader.
            variants = self._variants if self._variants else self._reader

            for record in variants:

                if only_snps and not record.is_snp:
                    continue

                if only_good and record.FILTER:
                    continue

                writer.write_record(record)
                written_variants += 1

        return written_variants

# --------------------------------------------------------------------------------------------------

    def write_to_json(self, vcf_file_name, verbose=False):
        """Write _variants to a json file.

        Parameters:
        -----------


        Returns:
        --------
        0
        """

        data = {'annotations': {}, 'positions': {}}

        for info_key, metadata in self.out_template.metadata.iteritems():
            data['annotations'][info_key] = metadata

        for record in self._variants:
            # Inject a property about uncallable genotypes into record.
            record.__setattr__("is_uncallable", is_uncallable(record))  # is_uncallable = types.MethodType(is_uncallable, record)

            # SKIP indels, if not handled then can cause REF base to be >1
            if record.is_indel and not (record.is_uncallable or record.is_monomorphic) or len(record.REF) > 1:
            # if len(record.REF) > 1:
                # print "%s\t%s\t%s\t%s\t%s" % (sample_name,record.POS,-1,record.FILTER,record)
                continue
                # if record.is_deletion and not record.is_uncallable:
                #    continue

            try:
                _ = data['positions'][record.CHROM]
            except KeyError:
                data['positions'][record.CHROM] = {"G": [], "A": [], "T": [], "C": [], "N": [], "-": []}

            # IF this is uncallable genotype, add gap "-"
            if record.is_uncallable:
                data['positions'][record.CHROM]["-"].append(record.POS)

            elif not record.FILTER:
                # If filter PASSED!
                if record.is_snp:
                    if len(record.ALT) > 1:
                        logging.info("POS %s passed filters but has multiple alleles REF: %s, ALT: %s. Inserting N", record.POS, str(record.REF), str(record.ALT))
                        data['positions'][record.CHROM]["N"].append(record.POS)

                    else:
                        data['positions'][record.CHROM][str(record.ALT[0]).upper()].append(record.POS)

            # Filter(s) failed
            elif record.is_snp:
                data['positions'][record.CHROM]["N"].append(record.POS)
            else:
                data['positions'][record.CHROM]["N"].append(record.POS)

        # print summary information if requested
        if verbose == True:
            logging.info("%s chromosomes were found in the VCF file" % len(data['positions']))
            for chromosome in data['positions']:
                logging.info("")
                logging.info("%s: " % chromosome)
                for position_type in ["G", "A", "T", "C", "N", "-"]:
                    logging.info("\t%s: %s" % (position_type, len(data['positions'][chromosome][position_type])))
                logging.info("-----------------")

        json_string = json.dumps(data)
        json_file_name = vcf_file_name.replace('.vcf', '.json.gz')

        with gzip.open(json_file_name, 'wb') as fJson:
            fJson.write(json_string)

        return 0

# --------------------------------------------------------------------------------------------------

    def _write_bad_variants(self, vcf_out):
        """**PRIVATE:** Write only those records that **haven't** passed."""
        written_variants = 0
        # Check if the output file ends with .gz, then compress data.
        open_func = gzip.open if vcf_out.endswith(".gz") else open
        with open_func(vcf_out, "w") as out_vcf:
            writer = vcf.Writer(out_vcf, self.out_template)
            for record in self._variants:
                if record.FILTER != "PASS" and record.FILTER is not None:
                    writer.write_record(record)
                    written_variants += 1
        return written_variants

    def _update_filters(self, new_filters):
        """Update internal filters in the output template."""
        for new_filter, filter_data in new_filters.items():
            self.out_template.filters[new_filter] = filter_data


class VariantCaller(PHEMetaData):
    """Abstract class used for access to the implemented variant callers."""

    __metaclass__ = abc.ABCMeta

    def __init__(self, cmd_options=None):
        """Constructor for variant caller.

        Parameters
        ----------
        cmd_options: str, optional
            Command options to pass to the variant command.
        """
        self.cmd_options = cmd_options

        super(VariantCaller, self).__init__()

    @abc.abstractmethod
    def make_vcf(self, *args, **kwargs):
        """Create a VCF from **BAM** file.

        Parameters
        ----------
        ref: str
            Path to the reference file.
        bam: str
            Path to the indexed **BAM** file for calling _variants.
        vcf_file: str
            path to the VCF file where data will be written to.
        make_aux: bool, optional
            True/False create auxilliary files (default: False).

        Returns
        -------
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

    def validate(self):
        """Check if the command can run."""
        if self.get_version() == "n/a":
            raise Exception("%s is not available in your PATH." % self.name)
