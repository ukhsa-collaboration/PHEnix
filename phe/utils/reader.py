'''
:Date: 3 June, 2016
:Author: Public Health England
'''
from collections import Counter, defaultdict

import vcf

from phe.utils import is_uncallable


class ParallelVCFReader(object):
    """Read multiple VCFs in parallel - one position at the time."""

    def __init__(self, vcfs):
        """Instantiate ParallelVCFReader.

        Parameters
        ----------
        vcfs: list
            List of path s to the VCF files to read.
        """
        self._readers = { vcf_in: vcf.Reader(filename=vcf_in) for vcf_in in vcfs}

        self._records = {}
        self.update()

    def __iter__(self):
        return self.get_records()

    def _vote_best_ref(self, refs):
        """Pick the best reference by counting most frequent occurence."""
        counts = 0
        try:
            counts = max(Counter(refs).iteritems(), key=lambda x: x[1])[0]
        except ValueError:
            counts = 0

        return counts

    def get_samples(self):
        """Get the list of sample names from the VCFs."""
        return [reader.samples[0] for reader in self._readers.itervalues()]

    def update(self, ids=None):
        """Update records in the readers. If *ids* is **not**
        specified, then all readers are updated. Otherwise, only
        those with the same *id* are updated.

        Parameters
        ----------
        ids: list, optional
            Optional list of IDs to update. If **None**, all will be updated.

        """
        if ids is None:
            ids = self._readers.keys()

        for k in ids:
            try:
                self._records[k] = self._readers[k].next()
                self._records[k].__setattr__("is_uncallable", is_uncallable(self._records[k]))
            except StopIteration:
                self._records[k] = None

    def get_records(self):
        """Generator of records, one position at the time. If a position has
        more than 1 record in any of the readers, a list of all records with that
        position is returned. The user can deal appropriately with the records list.

        Returns
        -------
        chrom: str
            Chromosome
        pos: int
            Position
        records: dict
            Records in dictionary of lists for each sample.
            
            .. code-block:: python

                {
                    "sample1": [record1, record2], 
                    "sample2": [record3]
                }

        """

        chrom = None
        # Iterate while we have reader and they are not None.
        while self._readers and not all(r is None for r in self._records.itervalues()):

            # Find the 'best' - most frequent chromosome.
            if chrom is None:
                chrom = self._vote_best_ref([ record.CHROM for record in self._records.itervalues() if record])

            # Find the minimum position for that chromosome.
            try:
                # Find the lowest record position in all readers.
                min_pos = min([ record.POS for record in self._records.itervalues() if record and record.CHROM == chrom])
            except ValueError:
                # This happens when all records have different chrom.
                chrom = self._vote_best_ref([ record.CHROM for record in self._records.itervalues() if record])
                continue

            # Only get the records for vcfs that are in lowest currect position.
            records = defaultdict(list)
            # Pull out all records for that position + chromosome from readers.
            for vcf_in, record in self._records.iteritems():

                while record and record.CHROM == chrom and record.POS == min_pos:
                    records[self._readers[vcf_in].samples[0]].append(record)
                    self.update([vcf_in])

                    record = self._records[vcf_in]

            if not records:
                continue

            # Yield them to iterator.
            yield chrom, min_pos, records

            # Update the records only of the readers we used, i.e. lowest position.
            for vcf_in in self._readers:
                if self._records[vcf_in] and self._records[vcf_in].POS == min_pos and self._records[vcf_in].CHROM == chrom:
                    self.update([vcf_in])
