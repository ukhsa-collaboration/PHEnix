==========
Annotators
==========

Individual VCFs can be annotated using custom annotators. Currently available annotators:

- **coverage** - Annotates with information about _mean_ and _dev_ of the depth in the VCF (using DP from INFO).

The data can be accessed from the metadata field in the VCF in the following way:

.. code:: python

   r = vcf.Reader(filename="/path/to/vcf")
   print r.metadata["coverageMetaData"][0]["mean"]
   print r.metadata["coverageMetaData"][0]["dev"]
