# PHEnix

Full documentation can be found at: http://phenix.readthedocs.io/en/latest/

## Installation

From source:

```bash
$ git clone https://github.com/phe-bioinformatics/PHEnix.git
$ pip2 install -e PHEnix
```

Directly from github:

```bash
$ pip install git+https://github.com/phe-bioinformatics/PHEnix.git
```

Installing from Pip:

```bash
$ pip install PHEnix
```

## Prereq

The reference needs to be index in an appropriate way for different mapper/variant callers. This can be done using:

```bash
$ phenix.py prepare_reference.py --mapper [bwa | bowtie2] --variant [gatk | mpileup] --reference <path_to_reference>
```

Internally, when the pipeline is ran internally, the indexing will be done automatically, because it
assumes that a reference is present for each sample. It has to be done differently for other
use cases because if multiple jobs are spawned and each one tries to index the same file, this may
lead to corrupted indices.

**N.B> The reference file can't have _fas_ extension, because Picard Tools throws an exception.**

## Calling SNPs

If you have a sample and you want to have one-stop-shop analysis run the following:

```bash
$ phenix.py run_snp_pipeline.py -r1 <path to R1.fastq> -r2 <path to R2.fastq> \
-r <path to reference> --sample-name <NAME> --mapper bwa --variant gatk \
--filters min_depth:5,mq_score:30
```

This will map with **BWA** and call variants with **GATK**. Intermediate files are written into the same directory you run this
command from. **--sample-name** option is very important, it specifies what output files will be called and the read group in the BAM
file. If you omit it, **test_sample** will be used.

## Filters

One of the key parts of the VCF processing is to filter quality calls. To do this we have created a flexible interface with
variety of filters available:

- **qual_score** - Filter records where **QUAL** score is below given threshold.
- **ad_ratio** - Filter, defined by **gatk**, records where ratio of alt allele to sum of all alleles is below given fraction.
- **dp4_ratio** - Similar to **ad_ratio**, but used in **mpileup** variant caller.
- **mq_score** - Filter records that fall below specified **MQ** score (from _INFO_ field).
- **mq0_ratio** - Filter, defined by **gatk**, records with **MQ0** to **DP** ratio _above_ given threshold (both values are from _INFO_ field).
- **mq0f_ratio** - Similar to the **mq0_ratio**, but used in **mpileup** variant caller.
- **qg_score** - Filter records that fall below specified **GQ** score (from **first** sample).
- **min_depth** - Filter records with mean depth below specified threshold (**DP** from sample field).
- **uncall_gt** - Filter records with uncallable genotypes in VCF (**GT** is ./.).

All filters are applied for each position and those positions that pass ALL filters are kept as quality calls. Positions
failing filter will be kept for future reference and creating fasta files, when needed. To specify filters to be used, simply
list them as key:threshold pairs separated by comma(,). For filters that don't require threshold, leave blank after ':'.

## Annotators

Individual VCFs can be annotated using custom annotators. Currently available annotators:

- **coverage** - Annotates with information about _mean_ and _dev_ of the depth in the VCF (using DP from INFO).

The data can be accessed from the metadata field in the VCF in the following way:

```python
r = vcf.Reader(filename="/path/to/vcf")
print r.metadata["coverageMetaData"][0]["mean"]
print r.metadata["coverageMetaData"][0]["dev"]
```

## Converting from VCF to FASTA

A lot of downstream applications take on FASTA formatted files, not VCFs. We have included a script for converting VCF data to
FASTA format.

```bash
$ phenix.py vcf2fasta -d <path to directory of VCFs> -o <path to output FASTA>
```

This tool is also able to filter out samples and columns that may be bad quality. E.g. **--sample-Ns** specifies maximum fraction of Ns
present in a sample. **--Ns** specifies maximum fraction of Ns allowed per column. **--with-mixtures** can specify if mixed position
should be output as over certain fraction. First, samples are sifted out then columns are checked. **--with-stats** allows to output
genomic positions of the called SNPs. Each line corresponds to a character position in the FASTA file.
E.g. to get the position of the 23rd SNP, go to 23rd line in the positions file. With this option, if numpy and matplotlib
are installed, **plots** directory will be created and a summary plot is generated, summarising called SNPs,
Ns, mixtures, and bases with 0 depth. Providing **--with-dist-matrix** will also write distance matrix to the specified path.

## Requirements

A lot of functionality depends on the presence of existing 3rd party tools:

Mappers:

- BWA
- Bowtie2

Variant Caller:

- GATK
- MPileup

In order for them to function properly, they need to be already in you *PATH*. For commands that
run through Java archives, please set appropriate environment variable (see below).

Python:

- Python >= 2.7
- argparse
- PyVCF
- PyYAML
- matplotlib (_optional_)
- bintrees (_optional_)
- numpy (_optional_)
- matplotlib.venn (_optional_)
- psycopg2 (_optional_)

### Cython

Both PyVCF and bintrees will install C bindings if Cython is install first.
These run faster than the native Python implementation. By installing it in
this order, you will avoid the following warnings from appearing when
running `phenix.py`:

```
Warning: FastBinaryTree not available, using Python version BinaryTree.
Warning: FastAVLTree not available, using Python version AVLTree.
Warning: FastRBTree not available, using Python version RBTree.
```

To avoid this, make sure that PyVCF and bintrees are not installed yet:

```
pip2 uninstall PyVCF
pip2 uninstall bintrees
pip2 install cython
pip2 install PyVCF bintrees
```

## 3rd Party Requirements

### Samtools

Samtools can be downloaded from https://github.com/samtools/samtools. It is used to filter and convert to SAM/BAM files and in mpileup variant caller.

### BCFTools

BCFtools can be downloaded from https://samtools.github.io/bcftools/. It is used for calling variants in mpileup.

### BWA

The BWA mapper can be downloaded from http://bio-bwa.sourceforge.net/.

### Bowtie2

The Bowtie2 mapper is available from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml.

### GATK
GATK is available from https://www.broadinstitute.org/gatk/. Please read the licencing information before using (https://www.broadinstitute.org/gatk/about/#licensing)

Set *GATK_JAR* - full path to the GATK Java archive.

### Picard
The Picard tool suite is available from http://broadinstitute.github.io/picard/.

Picard is needed for GATK to create dictionary of reference fasta.
Either set *PICARD_TOOLS_PATH* - path to directory where different Picard jars are (older versions) or set *PICARD_JAR* - path to **picard.jar**.
Older Picard distributions have many different jars (use first suggestion), where as newer versions have merged all into one jar file.

## Contributing to the project
If you want to add your favourite Mapper/VariantCaller or add another Filter please refer to the documentation: [here](http://phoenix.readthedocs.org/en/latest/api/modules.html)


## History
### Contributors
Aleksey Jironkin, Ali Al-Shahib, Anais Painset, Anthony Underwood, Georgia Kapatai, Philip Ashton, Rediat Tewolde, Richard Myers, Tim Dallman, Ulf Schaefer

### The (re)birth of PHEnix
In the beginning there were 2 pipelines to call SNPs: upstairs and downstairs.
One day we decided it was a bad idea to have 2 pipelines that did similar things.
The existing pipelines faded away, and a new pipeline was born revitalised and regenerated from the previous two, just like the ancient myth of the phoenix or phenix as it was known as in middle English (https://en.wikipedia.org/wiki/Phoenix_(mythology))
