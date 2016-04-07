'''Filter VCF on MQ filter.

:Date: 24 Sep, 2015
:Author: Alex Jironkin
'''

import argparse
import logging

from phe.variant_filters import PHEFilterBase


class MQ0FFilter(PHEFilterBase):
    '''Filter sites by MQ0F (Total Mapping Quality Zero Reads to DP) ratio.'''

    name = "MinMQ0F"
    _default_threshold = 0.05
    parameter = "mq0f_ratio"

    @classmethod
    def customize_parser(self, parser):
        arg_name = self.parameter.replace("_", "-")
        parser.add_argument("--%s" % arg_name, type=float, default=self._default_threshold,
                help="Filter sites below given MQ0F ratio (default: %s)" % self._default_threshold)

    def __init__(self, args):
        """Min Mapping Quality Zero constructor."""
        # This needs to happen first, because threshold is initialised here.
        super(MQ0FFilter, self).__init__(args)

        # Change the threshold to custom gq value.
        self.threshold = self._default_threshold
        if isinstance(args, argparse.Namespace):
            self.threshold = args.mq_score
        elif isinstance(args, dict):
            try:
                self.threshold = float(args.get(self.parameter))
            except (TypeError, ValueError):
                logging.error("Could not retrieve threshold from %s", args.get(self.parameter))
                logging.error("This parameter requires to be a float!")
                raise Exception("Could not create MQ0F filter from parameters: %s" % args)

    def __call__(self, record):
        """Filter a :py:class:`vcf.model._Record`."""

        good_record = self._check_record(record)

        if good_record is not True:
            return good_record
        record_mq = record.INFO.get("MQ0F")

        if record_mq is None or record_mq > self.threshold:
            # FIXME: when record_mq is None, i,e, error/missing, what do you do?
            return record_mq or False
        else:
            return None

    def short_desc(self):
        short_desc = self.__doc__ or ''

        if short_desc:
            short_desc = "%s (MQ0F > %s)" % (short_desc, self.threshold)

        return short_desc
