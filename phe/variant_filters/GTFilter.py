'''Filter VCF on GT filter parameter.

:Date: 24 Sep, 2015

:Author: Alex Jironkin
'''

import argparse
import logging

from phe.variant_filters import PHEFilterBase


class UncallableGTFilter(PHEFilterBase):
    '''Filter uncallable genotypes'''

    name = "UncallGT"
    _default_threshold = None
    parameter = "uncall_gt"

    @classmethod
    def customize_parser(self, parser):
        arg_name = self.parameter.replace("_", "-")
        parser.add_argument("--%s" % arg_name, type=str, default=self._default_threshold,
                help="Filter sites below given GQ score (default: %s)" % self._default_threshold)

    def __init__(self, args):
        """Min Depth constructor."""
        # This needs to happen first, because threshold is initialised here.
        super(UncallableGTFilter, self).__init__(args)

        # Change the threshold to custom gq value.
        self.threshold = self._default_threshold
        if isinstance(args, argparse.Namespace):
            self.threshold = args.gq_score
        elif isinstance(args, dict):
            try:
                self.threshold = str(args.get(self.parameter))
            except (TypeError, ValueError):
                logging.error("Could not retrieve threshold from %s", args.get(self.parameter))
                logging.error("This parameter requires to be a string!")
                self.threshold = None

    def __call__(self, record):
        """Filter a :py:class:`vcf.model._Record`."""

        if len(record.samples) > 1:
            logging.warn("More than 1 sample detected. Only first is considered.")

        try:
            record_gt = record.samples[0].data.GT
        except AttributeError:
            logging.warn("Could not retrieve GQ score POS %i", record.POS)
            record_gt = None

        if record_gt is None:
            return "./."
        else:
            return None

    def short_desc(self):
        short_desc = self.__doc__ or ''

        if short_desc:
            short_desc = "%s (GT != ./. )" % (short_desc)

        return short_desc

    def is_gap(self):
        return True

    def is_n(self):
        return False
