'''Filter VCF on depth of coverage.

:Date: 24 Sep, 2015
:Author: Alex Jironkin
'''
import argparse
import logging

from phe.variant_filters import PHEFilterBase


class DepthFilter(PHEFilterBase):
    """Filter sites by depth."""

    name = "MinDepth"
    _default_threshold = 5
    parameter = "min_depth"

    @classmethod
    def customize_parser(self, parser):
        arg_name = self.parameter.replace("_", "-")
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
            try:
                self.threshold = int(args.get(self.parameter))
            except (TypeError, ValueError):
                logging.error("Could not retrieve integer threshold from %s", args.get(self.parameter))
                raise Exception("Could not create depth filter from parameters: %s" % args)

    def __call__(self, record):
        """Filter a :py:class:`vcf.model._Record`."""

        if len(record.samples) > 1:
            logging.warn("Currently we only filter VCFs with 1 sample. Only first sample will be used.")
            logging.warn("This parameter requires to be an integer!")

        try:
            record_dp = record.samples[0].data.DP
        except AttributeError:
            logging.debug("Falling back to INFO DP")
            record_dp = record.INFO.get("DP")

        if record_dp is None or record_dp < self.threshold:
            return record_dp or False
        else:
            return None

    def short_desc(self):
        short_desc = self.__doc__ or ''

        if short_desc:
            short_desc = "%s (DP > %i)" % (short_desc, self.threshold)

        return short_desc
