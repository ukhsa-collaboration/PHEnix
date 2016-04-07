============
Calling SNPs
============

If you have a sample and you want to have one-stop-shop analysis run the
following: 

.. code-block:: bash 

   $ run_snp_pipeline.py \
   -r1 <path to R1.fastq> \
   -r2 <path to R2.fastq> \
   -r <path to reference> \
   --sample-name <NAME> \
   --mapper bwa \
   --variant gatk \
   --filters min_depth:5,mq_score:30 \
   -o <path to output directory>

This will map with **BWA** and call variants with **GATK**. Intermediate files are written into
the specified **output directory**. **--sample-name** option is very
important, it specifies what output files will be called and the read group in
the BAM file. If you omit it, **test_sample** will be used.

Config
------

If you have a pipeline, or want to run the same settings for different samples this process can be simplefied
by having a single config file. The file is in **YAML** format and takes common options that may be used for all samples:

- mapper - **string**
- mapper-options - **string**
- variant - **string**
- variant-options - **string**
- filters - **dictionary**
- annotators - **list**

Below is an example of config file:

.. code:: YAML

   mapper: bwa
   mapper-options: -t 4
   variant: gatk
   variant-options: --sample_ploidy 2 --genotype_likelihoods_model SNP -rf BadCigar -out_mode EMIT_ALL_SITES -nt 1
   filters:
     ad_ratio: 0.9
     min_depth: 5
     qual_score: 30
     mq_score: 30
   annotators:
     - coverage
