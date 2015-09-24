'''
Created on 24 Sep 2015

@author: alex
'''
import argparse

from phe.variant_filters import PHEFilterBase

class DepthFilter(PHEFilterBase):
    """Filter sites by depth."""

    name = "MinDepth"
    _default_threshold = 5
    _parameter = "min_depth"

    @classmethod
    def customize_parser(self, parser):
        arg_name = self._parameter.replace("_", "-")
        parser.add_argument("--" % arg_name, type=int, default=self._default_threshold,
                help="Filter sites below minimum depth (default: %s)" % self._default_threshold)

    def __init__(self, args):
        """Min Depth constructor."""
        # This needs to happen first, because threshold is initialised here.
        super(DepthFilter, self).__init__(args)

        # Change the threshold to custom dp value.
        self.threshold = self._default_threshold
        if isinstance(args, argparse.Namespace):
            self.threshold = args.min_depth
        elif isinstance(args, dict):
            self.threshold = args.get(self._parameter)

    def __call__(self, record):
        """Filter a :py:class:`vcf.model._Record`."""
        record_dp = record.INFO.get("DP")
        if record_dp is None or record_dp < self.threshold:
            return record_dp
        else:
            return None

    def short_desc(self):
        short_desc = self.__doc__ or ''

        if short_desc:
            short_desc = "%s (DP > %i)" % (short_desc, self.threshold)

        return short_desc
