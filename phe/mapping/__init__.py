"""

"""

import abc

class Mapper(object):
    """
    Blah
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, cmd_options=None):
        self._cmd_options = cmd_options

    @abc.abstractmethod
    def make_sam(self, **kwargs):
        raise NotImplementedError("make_same is not implemented yet.")

    @abc.abstractmethod
    def get_info(self):
        raise NotImplementedError("get_info is not implemented yet.")
