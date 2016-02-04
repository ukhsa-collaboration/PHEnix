
import re

import vcf
from vcf.model import _Record


class PHEVCFLightReader(vcf.Reader):

    def next(self):
        '''Return the next record in the file.'''
        line = self.reader.next()
        row = re.split(self._separator, line.rstrip())
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(self._parse_alt, row[4].split(','))

        try:
            qual = int(row[5])
        except ValueError:
            try:
                qual = float(row[5])
            except ValueError:
                qual = None

        filt = row[6]
        if filt == '.':
            filt = None
        elif filt == 'PASS':
            filt = []
        else:
            filt = filt.split(';')
#         info = self._parse_info(row[7])

        fmt = "GT:AD"
#         try:
#             fmt = row[8]
#         except IndexError:
#             fmt = None
#         else:
#             if fmt == '.':
#                 fmt = None

        record = _Record(chrom, pos, ID, ref, alt, qual, filt,
                None, fmt, self._sample_indexes)

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record
