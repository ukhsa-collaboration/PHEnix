'''
:Date: 10 Nov, 2015
:Author: Alex Jironkin
'''
import abc
from collections import OrderedDict
import glob
import inspect
import logging
import os
import sys

from phe.metadata import PHEMetaData

class Annotator(PHEMetaData):
    """Base class for Annotators"""
    __meta__ = abc.ABCMeta
    name = None

    def __init__(self, name):
        self.name = name

    @abc.abstractmethod
    def annotate(self, vcf_path=None):
        raise NotImplementedError("This is a base method to be implemented.")

    @abc.abstractmethod
    def get_meta_values(self):
        raise NotImplementedError("This is a base method for meta values that needs to be implented")

    def get_meta(self):
        od = self.get_meta_values()
        return OrderedDict({"%sMetaData" % self.name: [od]})

def dynamic_annotator_loader():
    """Fancy way of dynamically importing existing annotators.
    
    Returns
    -------
    dict:
        Available annotators dictionary. Keys are parameters that
        can be supplied to the filters.
    """

    # We assume the filters are in the same directory as THIS file.
    annotators_dir = os.path.dirname(__file__)
    annotators_dir = os.path.abspath(annotators_dir)

    # This is populated when the module is first imported.
    avail_annotators = {}

    # Add this directory to the syspath.
    sys.path.insert(0, annotators_dir)

    # Find all "py" files.
    for annotator_mod in glob.glob(os.path.join(annotators_dir, "*.py")):

        # Derive name of the module where filter is.
        filter_mod_file = os.path.basename(annotator_mod)

        # Ignore this file, obviously.
        if filter_mod_file.startswith("__init__"):
            continue

        # Import the module with a filter.
        mod = __import__(filter_mod_file.replace(".pyc", "").replace(".py", ""))

        # Find all the classes contained in this module.
        classes = inspect.getmembers(mod, inspect.isclass)
        for cls_name, cls in classes:
            # For each class, if it is a sublass of PHEFilterBase, add it.
            if cls_name != "Annotator" and issubclass(cls, Annotator):
                # The parameters are inherited and defined within each filter.
                avail_annotators[cls.name] = cls

    sys.path.remove(annotators_dir)

    return avail_annotators

_avail_annotators = dynamic_annotator_loader()

def available_annotators():
    """Return list of available filters."""
    return _avail_annotators.keys()

def make_annotators(config):
    """Create a list of annotators from *config*.
    
    Parameters
    ----------
    config: dict, optional
        Dictionary with name: value pairs. For each name, an
        appropriate Annotator will be found and instanciated.
        
    Returns
    -------
    list:
        List of :py:class:`Annotator` annotators.
    """
    annotators = []

    if config:
        for custom_annotator in config:
            if custom_annotator in _avail_annotators:
                annotators.append(_avail_annotators[custom_annotator]())
            else:
                logging.warn("Could not find appropriate filter for %s",
                             custom_annotator)

    return annotators

