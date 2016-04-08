"""
Classes and functions for working with variant filters.

.. _implementing-filter:

Implementing Filter
-------------------

The :py:class:`PHEFilterBase` class was designed to be implemented in order
to provide new filters for processing. As such it should be easy to extend
:py:class:`PHEFilterBase` to achieve this.

If you are not familiar with Python and classes please read
`this <http://www.tutorialspoint.com/python/python_classes_objects.htm>`_.
It is a general introduction to the subject.

Implementation of a filter closely follows the structure for
`PyVCF <http://pyvcf.readthedocs.org/en/latest/FILTERS.html#adding-a-filter>`_

PHEBaseFilter adds additional infomations that is also added to the VCF header:

- :py:attr:`PHEFilterBase.parameter` - Atrribute that allows dynmic loader to pickup the filter.This is the name used in command line arguments and configs.
- :py:attr:`PHEFilterBase._default_threshold` - Default threshold used for the filter when non specified.
- :py:meth:`PHEFilterBase.short_desc` - Method for defining what the short description is dynamically.
- :py:meth:`PHEFilterBase.__call__()` - Function that does the filtering (the same as PyVCF).

The most important of these is the **__call__** function, as this is what
does the filtering. It needs to return the following:

- **None** if the record **PASSES** the filter.
- Non **None** expressions are treated as failures.

.. NOTE::
   While it may not be obvious, but Python always returns a value.
   When no *return* statement is specified, **None** is returned
   to the calling function.

   **False** is returned when data can not be accessed by the filter.

Once you have implemented the filter and saved it in the *variant_filters*
directory, it should be automatically picked up by the system and available
to filter on. To verify this run:

.. code-block:: bash

    run_snp_pipeline.py --help

You filter should now be visible in the list of available filters.

Example
-------

.. code-block:: python

   class MyFilter(PHEFilterBase):
       name="CoolQUAL"
       _default_threshold=40
       parameter=cool_qual

       def __init__(self, args):
           # Construct any specific details for the object

       def __call__(self, record):
           if record.QUAL < self.threshold:
               return record.QUAL

       def short_desc(self):
           return "My cool QUAL filter"


:Date: 24 Sep, 2015
:Author: Alex Jironkin
"""

from __builtin__ import __import__
from abc import abstractproperty
import abc
import argparse
import glob
import inspect
import logging
import os
import re
import sys

import vcf
import vcf.filters
from vcf.parser import _Filter

IUPAC_CODES = {frozenset(["A"]): "A",
               frozenset(["C"]): "C",
               frozenset(["G"]): "G",
               frozenset(["T"]): "T",
               frozenset(["A", "G"]): "R",
                frozenset(["C", "T"]): "Y",
                frozenset(["G", "C"]): "S",
                frozenset(["A", "T"]): "W",
                frozenset(["G", "T"]): "K",
                frozenset(["A", "C"]): "M",
                frozenset(["C", "G", "T"]): "B",
                frozenset(["A", "G", "T"]): "D",
                frozenset(["A", "C", "T"]): "H",
                frozenset(["A", "C", "G"]): "V"
              }

class PHEFilterBase(vcf.filters.Base):
    """Base class for VCF filters."""
    __meta__ = abc.ABCMeta

    magic_sep = ":"
    decoder_pattern = re.compile(magic_sep)

    @abc.abstractproperty
    def parameter(self):
        """Short name of parameter being filtered."""
        return self.parameter

    @abc.abstractproperty
    def _default_threshold(self):
        """Default threshold for filtering."""
        return self._default_threshold

    def __init__(self, args):
        super(PHEFilterBase, self).__init__(args)

        # Change the threshold to custom gq value.
        self.threshold = self._default_threshold

        if isinstance(args, dict):
            self.threshold = args.get(self.parameter)

    def __str__(self):
        return self.filter_name()

    def _check_record(self, record):
        if self.is_uncallable(record):
            return False
        elif record.is_monomorphic:
            return None
        else:
            return True

    @abc.abstractmethod
    def short_desc(self):
        """Short description of the filter (included in VCF)."""
        raise NotImplementedError("Get short description is not implemented.")

    def get_config(self):
        """This is used for reconstructing filter."""
        return {self.parameter: self.threshold}

    def filter_name(self):
        """Create filter names by their parameter separated by magic.
        E.g. if filter parameter is ad_ratio and threshold is 0.9 then
        ad_ratio:0.9 if the filter name.
        """
        return "%s%s%s" % (self.parameter, self.magic_sep, self.threshold)

    @staticmethod
    def decode(filter_id):
        """Decode name of filter."""
        conf = {}

        if PHEFilterBase.magic_sep in filter_id:
            info = PHEFilterBase.decoder_pattern.split(filter_id)
            assert len(info) == 2
            conf[info[0]] = info[1]
        return conf

    def is_gap(self):
        return False

    def is_n(self):
        return True

    def is_uncallable(self, record):
        """Filter a :py:class:`vcf.model._Record`."""

        if len(record.samples) > 1:
            logging.warn("More than 1 sample detected. Only first is considered.")

        try:
            record_gt = record.samples[0].data.GT
        except AttributeError:
            logging.warn("Could not retrieve GQ score POS %i", record.POS)
            record_gt = "./."

        if record_gt == "./.":
            return True
        else:
            return False

    @staticmethod
    def call_concensus(record):
        if not record.FILTER:
            sample_ad = frozenset([str(c).upper() for c in record.ALT])
            return IUPAC_CODES.get(sample_ad, "N")

        else:
            sample_ad = frozenset([str(c).upper() for c in record.ALT] + [record.REF])

            return IUPAC_CODES.get(sample_ad, "N")

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
                avail_filters[cls.parameter] = cls

    sys.path.remove(filter_dir)

    return avail_filters

_avail_filters = dynamic_filter_loader()

def available_filters():
    """Return list of available filters."""
    return _avail_filters.keys()

def str_to_filters(filters):
    """Convert from filter string to array of filters.
    E.g. ad_ration:0.9,min_depth:5
    
    Parameters
    ----------
    filters: str
        String version of filters, separated by comma.
    
    Returns
    -------
    list:
        List of :py:class:`phe.variant_filters.PHEFilterBase` instances.
    """

    config = {}
    for kv_pair in filters.split(","):
        pair = kv_pair.split(":")
        assert len(pair) == 2, "Filters should be separated by ':' %s" % kv_pair

        # We don't care about casting them to correct type because Filters
        #    will do it for us.
        config[pair[0]] = pair[1]

    return make_filters(config)

def make_filters(config):
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
    filters = []

    if config:
        for custom_filter in config:
            if custom_filter in _avail_filters:
                filters.append(_avail_filters[custom_filter](config))
            else:
                logging.error("Could not find appropriate filter for %s",
                             custom_filter)
                raise Exception("Filter %s could not be created. Please check your filter arguments." % custom_filter)

    return filters
