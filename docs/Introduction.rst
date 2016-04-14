============
Introduction
============
:Author: Public Health England

:Date: |today|

Installation
------------

From source:

.. code-block:: bash

   git clone https://github.com/phe-bioinformatics/PHEnix.git
   pip2 install -e PHEnix
   
Directly from github:

.. code-block:: bash
   
   pip install git+git://github.com/phe-bioinformatics/PHEnix.git
   
.. NOTE:: Installing from Pip - Coming Soon.


Overview
--------

This code was designed to allow users to input fastq files and a reference sequence and perform:

 - Reference mapping
 - VCF generation
 - VCF filtering
 - FASTA sequence of SNPs
 
The process is comprised of three steps:

 - Reference sequence preparation (prepare_reference.py)
 - Mapping and filtered VCF generation (run_snp_pipeline.py)
 - FASTA generation from single or multiple VCFs (vcf2fasta.py)

**Example:**

prepare a bwa and gatk reference using a fasta file (myref.fasta)

.. code-block:: bash

   prepare_reference.py \
   --mapper bwa \
   --variant gatk \
   --reference myref.fasta
   
map, call and filter variants on fastq files (my.R1.fastq, my.R2.fastq). Filter SNPs on minimum depth, mapping quality and AD ratio

.. code-block:: bash

   run_snp_pipeline.py \
   -r1 my.R1.fastq \
   -r2 my.R2.fastq \
   -r myref.fasta \
   --sample-name mysample \
   --mapper bwa \
   --variant gatk \
   --filters min_depth:5,mq_score:30,ad_ratio:0.9
   
generate a FASTA file of SNPs using filtered VCFs in current directory

.. code-block:: bash

   vcf2fastq.py -d ./ -o output.fasta --regex filtered

Requirements
------------

A lot of functionality depends on the presence of existing 3rd party tools:

Mappers:

* BWA - Download from [https://github.com/lh3/bwa]

* Bowtie2 - Download from [https://github.com/BenLangmead/bowtie2]



Variant Caller:

- GATK - Download from [https://www.broadinstitute.org/gatk/download/]
  * Picard - Download from [http://broadinstitute.github.io/picard/]

- MPileup - Download from [https://github.com/samtools/samtools]

In order for them to function properly, they need to be already in you **PATH**. For commands that run through Java archives, please set appropriate environment variable (see below).

Python
------

- Python >= 2.7

- argparse

- PyVCF

- PyYAML

- matplotlib (_optional_)

- bintrees (_optional_) - If you are using **vcf2fasta.py**

- numpy (_optional_)

- matplotlib.venn (_optional_) - psycopg2 (_optional_)


3rd Party Requirements
----------------------



Samtools
--------

Samtools Samtools can be downloaded from https://github.com/samtools/samtools. It is used to filter and convert to SAM/BAM files and in mpileup variant caller.


BCFTools
---------

BCFtools can be downloaded from https://github.com/samtools/bcftools. It is used for calling variants in mpileup. 

BWA Heng Li's mapper can be downloaded from https://github.com/lh3/bwa. 

Bowtie2 Bowtie2 mapper available from https://github.com/BenLangmead/bowtie2. 

GATK Set *GATK_JAR* - full path to the GATK Java archive.

Picard Tools
------------

Picard is needed for GATK to create dictionary of reference fasta. Either set *PICARD_TOOLS_PATH* - path to directory where different Picard jars are or set *PICARD_JAR* - path to **picard.jar**. Older Picard distributions have many different jars (use first suggestion), where as newer versions have merged all into one jar file. 


