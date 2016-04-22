'''
:Date: 22 Sep, 2015

:Author: Alex Jironkin
'''
import glob
import inspect
import logging
import os
import sys

from phe.mapping import Mapper

def dynamic_mapper_loader():
    """Fancy way of dynamically importing existing mappers.
    
    Returns
    -------
    dict:
        Available mappers dictionary. Keys are names that
        can be used for existing mapper implementations.
    """

    # We assume the mappers are in the same directory as THIS file.
    mappers_dir = os.path.dirname(__file__)
    mappers_dir = os.path.abspath(mappers_dir)

    # This is populated when the module is first imported.
    avail_mappers = {}

    # Add this directory to the syspath.
    sys.path.insert(0, mappers_dir)

    # Find all "py" files.
    for mapper_mod in glob.glob(os.path.join(mappers_dir, "*.py")):

        # Derive name of the module where mapper is.
        mapper_mod_file = os.path.basename(mapper_mod)

        # Ignore this __init__, only base class is there.
        if mapper_mod_file.startswith("__init__"):
            continue

        # Import the module with a mapper.
        mod = __import__(mapper_mod_file.replace(".pyc", "").replace(".py", ""))

        # Find all the classes contained in this module.
        classes = inspect.getmembers(mod, inspect.isclass)
        for cls_name, cls in classes:
            # For each class, if it is a sublass of Mapper, add it.
            if cls_name != "Mapper" and issubclass(cls, Mapper):
                # The name is inherited and defined within each mapper.
                avail_mappers[cls.name] = cls

    sys.path.remove(mappers_dir)

    return avail_mappers

_avail_mappers = dynamic_mapper_loader()

def available_mappers():
    """Return a list of available mappers."""
    return _avail_mappers.keys()

def factory(mapper=None, custom_options=None):
    """Create an instance of a mapper from _config_.
    
    Parameters
    ----------
    mapper: str, optional
        Name of the mapper to initialise.
    custom_options: str, optional
        Custom options to be passed to the mapper.
        
    Returns
    -------
    :py:class:`phe.mapping.Mapper`:
        Instance of a :py:class:`phe.mapping.Mapper` with provided name, or
        None if the mapper could not be found.
    """
    if mapper is not None and isinstance(mapper, str):

        mapper = mapper.lower()
        if mapper in _avail_mappers:
            try:
                obj = _avail_mappers[mapper](cmd_options=custom_options)
            except Exception as e:
                logging.error("Could not initialise mapper: %s", mapper)
                logging.error(e.message)
                obj = None

        else:
            logging.error("No implementation for %s mapper.", mapper)
            obj = None
        return obj

    logging.warn("Unknown parameters. Mapper could not be initialised.")
    return None
