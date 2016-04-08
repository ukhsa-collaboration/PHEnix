'''Filter VCFs on AD ratio.

:Date: 24 Sep, 2015
:Author: Alex Jironkin
'''

import argparse
import logging

from phe.variant_filters import PHEFilterBase


class ADFilter(PHEFilterBase):
    '''Filter sites by AD ratio.'''


    name = "ADRatio"
    _default_threshold = 0.9
    parameter = "ad_ratio"

    @classmethod
    def customize_parser(self, parser):
        arg_name = self.parameter.replace("_", "-")
        parser.add_argument("--%s" % arg_name, type=float, default=self._default_threshold,
                help="Filter sites below minimum ad ratio (default: %s)" % self._default_threshold)

    def __init__(self, args):
        """AD Ratio constructor."""
        # Call parent's constructor first.
        super(ADFilter, self).__init__(args)

        # Change the threshold to custom dp value.
        self.threshold = self._default_threshold
        if isinstance(args, argparse.Namespace):
            self.threshold = args.ad_ratio
        elif isinstance(args, dict):
            try:
                self.threshold = float(args.get(self.parameter))
            except (TypeError, ValueError):
                logging.error("Could not retrieve threshold from %s", args.get(self.parameter))
                logging.error("This parameter requires to be a float!")
                raise Exception("Could not create AD filter from parameters: %s" % args)


    def __call__(self, record):
        """Filter a :py:class:`vcf.model._Record`."""

        good_record = self._check_record(record)

        if good_record is not True:
            return good_record

        if len(record.samples) > 1:
            logging.warn("More than 1 sample detected. Only first is considered.")

        try:
            record_ad = record.samples[0].data.AD

            # FIXME: when record length is > 2, what do you do?
            assert len(record_ad) == 2, "AD data is incomplete POS: %i" % record.POS

            depth = sum(record.samples[0].data.AD)

            ratio = float(record_ad[1]) / depth
        except Exception:
            logging.warn("Could not calculate AD ratio from %s POS: %s", record, record.POS)
            ratio = None

        if ratio is None or ratio < self.threshold:
            # FIXME: When ratio is None, i.e. error, what do you do?
            return ratio or False
        else:
            return None

    def short_desc(self):
        short_desc = self.__doc__ or ''

        if short_desc:
            short_desc = "%s (AD ratio > %s )" % (short_desc, self.threshold)

        return short_desc
