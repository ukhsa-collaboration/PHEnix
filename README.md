# Phoenix

If you have a sample and you want to have one-stop-shop analysis run the following:

```run_snp_pipeline.py -r1 <path to R1.fastq> -r2 <path to R2.fastq> -r <path to reference> --sample-name --mapper bwa --variant gatk --filters min_depth:5,mq_score:30```

This will map with **BWA** and call variants with **GATK**. Intermediate files are written into the same directory you run this 
command from. **--sample-name** option is very important, it specifies what output files will be called and the read group in the BAM
file. If you omit it, **test_sample** will be used.

## Filters

One of the key parts of the VCF processing is to filter quality calls. To do this we have created a flexible interface with
variaety of filters avaialble:

- **qual_score** - Filter records where QUAL score is below given threshold.
- **ad_ration** - Filter records where allele ration is below given fraction.
- **mq_score** - Filter records that fall below specified MQ score.
- **mq0_score** - Filter records with MQ0 to DP ratio _above_ given threshold.
- **qg_score** - Filter records that fall below specified GQ score.
- **min_depth** - Filter records with mean depth below specified threshold.
- **uncall_gt** - Filter records with uncallable genotypes in VCF.

All filters are applied for each position and those positions that pass ALL filters are kept as quality calls. Positions
failing filter will be kept for future reference and creating fasta files, when needed. To specify filters to be used, simply
list them as key:threshold pairs separated by comma(,). For filters that don't require threshold, leave blank after ':'. 

## Requirements

A lot of functionality depends on the presence of existing 3rd party tools:

Mappers:

- BWA 

Variant Caller:

- GATK

In order for them to function properly, they need to be already in you *PATH*. For commands that
run through Java archives, please set appropriate environment variable (see below).

Python:

- argparse
- PyVCF
- PyYAML

## 3rd Party Requirements

### Samtools

Samtools can be downloadded from blurb. It is used to filter and convert to SAM/BAM files.

### BWA

Heng Li's mapper can be downloaded from blurb.

### GATK
Set *GATK_JAR* - full path to the GATK Java archive.

### Picard
Picard is needed for GATK to create dictionary of reference fasta.
Either set *PICARD_TOOLS_PATH* - path to directory where different Picard jars are or set *PICARD_JAR* - path to **picard.jar**.
Older Picard distributions have many different jars (use first suggestion), where as newer versions have merged all into one jar file. 



## History

In the beginning there were 2 pipelines to call SNPs: upstairs and downstairs.
One day we decided it was a bad idea to have 2 pieplines that did similar things.
The existing pipelines died quiet software death and new pipeline was born, to rule
them all.