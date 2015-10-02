'''
Created on 22 Sep 2015

@author: alex
'''
import glob
import inspect
import logging
import os
import sys

from phe.variant import VariantCaller

def dynamic_caller_loader():
    """Fancy way of dynamically importing existing variant callers.
    
    Returns
    -------
    dict:
        Available variant callers dictionary. Keys are parameters that
        can be used to call variants.
    """

    # We assume the caller are in the same directory as THIS file.
    variants_dir = os.path.dirname(__file__)
    variants_dir = os.path.abspath(variants_dir)

    # This is populated when the module is first imported.
    avail_callers = {}

    # Add this directory to the syspath.
    sys.path.insert(0, variants_dir)

    # Find all "py" files.
    for caller_mod in glob.glob(os.path.join(variants_dir, "*.py")):

        # Derive name of the module where caller is.
        caller_mod_file = os.path.basename(caller_mod)

        # Ignore __init__ file, only base class is there.
        if caller_mod_file.startswith("__init__"):
            continue

        # Import the module with a caller.
        mod = __import__(caller_mod_file.replace(".pyc", "").replace(".py", ""))

        # Find all the classes contained in this module.
        classes = inspect.getmembers(mod, inspect.isclass)
        for cls_name, cls in classes:
            # For each class, if it is a sublass of VariantCaller, add it.
            if cls_name != "VariantCaller" and issubclass(cls, VariantCaller):
                # The name is inherited and defined within each caller.
                avail_callers[cls.name] = cls

    sys.path.remove(variants_dir)

    return avail_callers

_avail_variant_callers = dynamic_caller_loader()

def available_callers():
    return _avail_variant_callers.keys()

def factory(config=None, variant=None, custom_options=None):
    if variant is not None and isinstance(variant, str):

        variant = variant.lower()
        if variant in _avail_variant_callers:
            return _avail_variant_callers[variant](cmd_options=custom_options)
        else:
            logging.error("No implementation for %s mapper.")
            return None

    logging.warn("Unknown parameters. Mapper could not be initialised.")
    return None
