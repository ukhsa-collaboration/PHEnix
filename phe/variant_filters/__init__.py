from __builtin__ import __import__
from abc import abstractproperty
import abc
import argparse
import glob
import inspect
import logging
import os
import sys

import vcf
import vcf.filters
from vcf.parser import _Filter

class PHEFilterBase(vcf.filters.Base):

    __meta__ = abc.ABCMeta

    @abc.abstractproperty
    def _parameter(self):
        return self._parameter

    @abc.abstractproperty
    def _default_threshold(self):
        return self._default_threshold

    def __init__(self, args):
        super(PHEFilterBase, self).__init__(args)

        # Change the threshold to custom gq value.
        self.threshold = self._default_threshold

        if isinstance(args, dict):
            self.threshold = args.get(self._parameter)

    @abc.abstractmethod
    def short_desc(self):
        raise NotImplementedError("Get short description is not implemented.")


def make_filters(*args, **kwargs):
    """Create a list of filters from *config*.
    
    Parameters:
    -----------
    config: dict, optional
        Dictionary with parameter: value pairs. For each parameter, an
        appropriate Filter will be found and instanciated.
        
    Returns:
    --------
    list:
        List of :py:class:`PHEFilterBase` filters.
    """
    config = kwargs.get("config")
    avail_filters = dynamic_filter_loader()

    filters = []

    if config:
        for custom_filter in config:
            if custom_filter in avail_filters:
                filters.append(avail_filters[custom_filter](config))
            else:
                logging.warn("Could not find appropriate filter for %s",
                             custom_filter)

    return filters

def dynamic_filter_loader():
    """Fancy way of dynamically importing existing filters.
    
    Returns
    -------
    dict:
        Available filters dictionary. Keys are parameters that
        can be supplied to the filters.
    """

    # We assume the filters are in the same directory as THIS file.
    filter_dir = os.path.dirname(__file__)
    filter_dir = os.path.abspath(filter_dir)

    # This is populated when the module is first imported.
    avail_filters = {}

    # Add this directory to the syspath.
    sys.path.insert(0, filter_dir)

    # Find all "py" files.
    for filter_mod in glob.glob(os.path.join(filter_dir, "*.py")):

        # Derive name of the module where filter is.
        filter_mod_file = os.path.basename(filter_mod)

        # Ignore this file, obviously.
        if filter_mod_file.startswith("__init__"):
            continue

        # Import the module with a filter.
        mod = __import__(filter_mod_file.replace(".pyc", "").replace(".py", ""))

        # Find all the classes contained in this module.
        classes = inspect.getmembers(mod, inspect.isclass)
        for cls_name, cls in classes:
            # For each class, if it is a sublass of PHEFilterBase, add it.
            if cls_name != "PHEFilterBase" and issubclass(cls, PHEFilterBase):
                # The parameters are inherited and defined within each filter.
                avail_filters[cls._parameter] = cls

    sys.path.remove(filter_dir)

    return avail_filters

def filter_vcf(vcf_in, filters):
    """WIP This should go somewhere else."""

    # Make filters.
    good_records = []
    bad_records = []

    # Create a reader class from input VCF.
    reader = vcf.Reader(filename=vcf_in)

    # Add each filter we are going to use to the record.
    # This is needed for writing out proper #FILTER header in VCF.
    for record_filter in filters:
        # We know that each filter has short description method.
        short_doc = record_filter.short_desc()
        short_doc = short_doc.split('\n')[0].lstrip()

        reader.filters[record_filter.filter_name()] = _Filter(record_filter.filter_name(), short_doc)

    # For each record (POSITION) apply set of filters.
    for record in reader:
        for record_filter in filters:

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
            good_records.append(record)
        else:
            bad_records.append(record)

    return good_records, bad_records, reader
