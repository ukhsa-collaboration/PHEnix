'''
Created on 22 Sep 2015

@author: alex
'''
import logging

from phe.variant.GATKVariantCaller import GATKVariantCaller


_variant_map = {"gatk": GATKVariantCaller}

def factory(config=None, variant=None, threahds=1):
    if variant is not None and isinstance(variant, str):

        variant = variant.lower()
        if variant in _variant_map:
            return _variant_map[variant]()
        else:
            logging.error("No implementation for %s mapper.")
            return None

    logging.warn("Unknown parameters. Mapper could not be initialised.")
    return None