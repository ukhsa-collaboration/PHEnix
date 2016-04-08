=============
Prerequisites
=============

The refence needs to be index in an appropriate way for different mapper/variant callers. This can be done manually for specific tools or using:

.. code:: bash 

   prepare_reference.py --mapper [bwa | bowtie2] \
   --variant [gatk | mpileup] \
   --reference <path_to_reference> 
   
Internally, for the pipeline, the indexing will be done automatically, because software assumes that a reference is present for each sample. It has to be done differently for other use cases because if multiple processes trying to index the same file, this may lead to corrupted indeces.

.. NOTE:: The reference file can't have *fas* extension, because Picard Tools throws an exception.