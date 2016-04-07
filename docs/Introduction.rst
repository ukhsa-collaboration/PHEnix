============
Introduction
============
:Author: Public Health England

:Date: |today|


Requirements
------------

A lot of functionality depends on the presence of existing 3rd party tools:

Mappers:

* BWA - Download from [blob]

* Bowtie2 - Download from [blob]



Variant Caller:

- GATK - Download from []
  * Picard - Download from

- MPileup - Download from []

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

Samtools Samtools can be downloadded from blurb. It is used to filter and convert to SAM/BAM files and in mpileup variant caller.


BCFTools
---------

BCFtools can be downloaded from blurb. It is used for calling variants in mpileup. 

BWA Heng Li's mapper can be downloaded from blurb. 

Bowtie2 Bowtie2 mapper available from blurb. 

GATK Set *GATK_JAR* - full path to the GATK Java archive.

Picard Tools
------------

Picard is needed for GATK to create dictionary of reference fasta. Either set *PICARD_TOOLS_PATH* - path to directory where different Picard jars are or set *PICARD_JAR* - path to **picard.jar**. Older Picard distributions have many different jars (use first suggestion), where as newer versions have merged all into one jar file. 


