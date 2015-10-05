"""Classes and functions for working with variant filters."""

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


class PHEFilterBase(vcf.filters.Base):
    """Base class for VCF filters."""
    __meta__ = abc.ABCMeta

    magic_sep = ":"

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
            info = filter_id.split(PHEFilterBase.magic_sep)
            assert len(info) == 2
            conf[info[0]] = info[1]
        return conf

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

def make_filters(**kwargs):
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

    filters = []

    if config:
        for custom_filter in config:
            if custom_filter in _avail_filters:
                filters.append(_avail_filters[custom_filter](config))
            else:
                logging.warn("Could not find appropriate filter for %s",
                             custom_filter)

    return filters
