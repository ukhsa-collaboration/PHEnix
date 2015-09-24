'''
Created on 24 Sep 2015

@author: alex
'''

import argparse
import logging

from phe.variant_filters import PHEFilterBase


class GQFilter(PHEFilterBase):
    '''Filter sites by GQ score.'''

    name = "MinGQ"
    _default_threshold = 0
    _parameter = "gq_score"

    @classmethod
    def customize_parser(self, parser):
        arg_name = self._parameter.replace("_", "-")
        parser.add_argument("--%s" % arg_name, type=int, default=self._default_threshold,
                help="Filter sites below given GQ score (default: %s)" % self._default_threshold)

    def __init__(self, args):
        """Min Depth constructor."""
        # This needs to happen first, because threshold is initialised here.
        super(GQFilter, self).__init__(args)

        # Change the threshold to custom gq value.
        self.threshold = self._default_threshold
        if isinstance(args, argparse.Namespace):
            self.threshold = args.gq_score
        elif isinstance(args, dict):
            self.threshold = args.get(self._parameter)

    def __call__(self, record):
        """Filter a :py:class:`vcf.model._Record`."""

        if len(record.samples) > 1:
            logging.warn("More than 1 sample detected. Only first is considered.")

        try:
            record_gq = record.samples[0].data.GQ
        except Exception:
            logging.error("Could not retrieve GQ score POS %i", record.POS)
            record_gq = None

        if record_gq is None or record_gq < self.threshold:
            # FIXME: when record_gq is None, i,e, error, what do you do?
            return record_gq or ""
        else:
            return None

    def short_desc(self):
        short_desc = self.__doc__ or ''

        if short_desc:
            short_desc = "%s (GQ > %s)" % (short_desc, self.threshold)

        return short_desc
