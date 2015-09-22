'''
Created on 22 Sep 2015

@author: alex
'''
import logging

from phe.mapping.BWAMapper import BWAMapper


_mapper_map = {"bwa": BWAMapper}

def factory(config=None, mapper=None, threahds=1):
    if mapper is not None and isinstance(mapper, str):

        mapper = mapper.lower()
        if mapper in _mapper_map:
            return _mapper_map[mapper]()
        else:
            logging.error("No implementation for %s mapper.")
            return None

    logging.warn("Unknown parameters. Mapper could not be initialised.")
    return None