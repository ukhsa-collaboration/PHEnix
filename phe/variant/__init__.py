import abc


class VariantCaller(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        return self.name

    def __init__(self, cmd_options=None):
        pass

    @abc.abstractmethod
    def make_vcf(self, *args, **kwargs):
        raise NotImplementedError("make_vcf is not implemented yet.")
