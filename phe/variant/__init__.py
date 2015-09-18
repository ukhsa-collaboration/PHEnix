import abc


class VariantCaller(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, cmd_options=None):
        pass

    @abc.abstractmethod
    def make_vcf(self):
        raise NotImplementedError("make_vcf is not implemented yet.")
