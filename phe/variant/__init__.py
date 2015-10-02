import abc
import pickle

from vcf import filters
import vcf
from vcf.parser import _Filter

from phe.variant_filters import make_filters

class VCFTemplate(object):

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

    def __init__(self, vcf_in, filters=None):

        self.vcf_in = vcf_in

        self.out_template = None

        self.filters = make_filters(config=filters)
        self.good_vars = []
        self.bad_vars = []

    def make_variant_set(self, keep_only_snps=True):
        """WIP This should go somewhere else."""

        # Create a reader class from input VCF.
        reader = vcf.Reader(filename=self.vcf_in)

        # Add each filter we are going to use to the record.
        # This is needed for writing out proper #FILTER header in VCF.
        for record_filter in self.filters:
            # We know that each filter has short description method.
            short_doc = record_filter.short_desc()
            short_doc = short_doc.split('\n')[0].lstrip()

            reader.filters[record_filter.filter_name()] = _Filter(record_filter.filter_name(), short_doc)

        # For each record (POSITION) apply set of filters.
        for record in reader:
            for record_filter in self.filters:

                # Call to __call__ method in each filter.
                result = record_filter(record)

                # Record is KEPT if filter returns None
                if result == None:
                    continue

                # save some work by skipping the rest of the code
    #             if drop_filtered:
    #                 output_record = False
    #                 break

                # If we got this far, then record is filtered OUT.
                record.add_filter(record_filter.filter_name())
    #             if short_circuit: break

    #         if output_record:
    #             # use PASS only if other filter names appear in the FILTER column
    #             # FIXME: is this good idea?
    #             if record.FILTER is None or record.FILTER = 'PASS':
    #             output.write_record(record)
            # After applying all filters, check if FILTER is None.
            # If it is, then record PASSED all filters.
            if record.FILTER is None:
                record.FILTER = 'PASS'
                # FIXME: Does this work for indels?
                if keep_only_snps and record.is_snp:
                    self.good_vars.append(record)
            else:
                self.bad_vars.append(record)

        self.out_template = VCFTemplate(reader)

        return self.good_vars


    def write_variants(self, vcf_out, only_snps=False):
        # Write records to VCF file.
#         with open("filtered_bad.vcf", "w") as out_vcf:
#             writer = vcf.Writer(out_vcf, template)
#             for record in bad_var:
#                 writer.write_record(record)

        with open(vcf_out, "w") as out_vcf:
            writer = vcf.Writer(out_vcf, self.out_template)
            for record in self.good_vars:

                if only_snps and not record.is_snp:
                    continue

                writer.write_record(record)

    def _write_bad_variants(self, vcf_out):
        with open(vcf_out, "w") as out_vcf:
            writer = vcf.Writer(out_vcf, self.out_template)
            for record in self.bad_vars:
                writer.write_record(record)

    def serialise(self, out_file):

        filter_conf = {}
        for custom_filter in self.filters:
            filter_conf.update(custom_filter.get_config())

        d = {"template": self.out_template,
             "filter_conf": filter_conf,
             "good_vars": self.good_vars,
             "bad_vars": self.bad_vars
             }

        with open(out_file, "wb") as fp:
            pickle.dump(d, fp, -1)

class VariantCaller(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, cmd_options=None):
        self.cmd_options = cmd_options

    @abc.abstractmethod
    def make_vcf(self, *args, **kwargs):
        raise NotImplementedError("make_vcf is not implemented yet.")
