"""

"""

import abc

class Mapper(object):
    """
    Blah
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        return self.name

    def __init__(self, cmd_options=None):
        self.cmd_options = cmd_options

    @abc.abstractmethod
    def make_sam(self, *args, **kwargs):
        raise NotImplementedError("make_sam is not implemented yet.")

    @abc.abstractmethod
    def make_bam(self, *args, **kwargs):
        raise NotImplementedError("make_bam is not implemented yet.")

    @abc.abstractmethod
    def get_info(self):
        raise NotImplementedError("get_info is not implemented yet.")
