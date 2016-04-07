"""Metadata related information.

:Date: 22 Sep, 2015
:Author: Alex Jironkin
"""

import abc

class PHEMetaData(object):
    """Abstract class to provide interface for meta-data creation."""

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass

    @abc.abstractmethod
    def get_meta(self):
        """Get the metadata."""
        raise NotImplementedError("get meta has not been implemented yet.")
