=======
Filters
=======

One of the key parts of the VCF processing is to filter quality calls. To do this we have created a flexible interface with variaety of filters avaialble:


- **qual_score** - Filter records where **QUAL** score is below given threshold.

- **ad_ratio** - Filter, defined by **gatk**, records where ratio of alt allele to sum of all alleles is below given fraction.

- **dp4_ratio** - Similar to **ad_ratio**, but used in **mpileup** variant caller.

- **mq_score**- Filter records that fall below specified **MQ** score (from _INFO_ field).

- **mq0_ratio** - Filter, defined by **gatk**, records with **MQ0** to **DP** ratio _above_ given threshold (both values are from _INFO_ field).

- **mq0f_ratio** - Similar to the **mq0_ratio**, but used in **mpileup** variant caller.

- **qg_score** - Filter records that fall below specified **GQ** score (from **first** sample).

- **min_depth** - Filter records with mean depth below specified threshold (**DP** from sample field).

- **uncall_gt** - Filter records with uncallable genotypes in VCF (**GT** is ./.).

All filters are applied for each position and those positions that pass ALL filters are kept as quality calls. Positions failing filter will be kept for future reference and creating fasta files, when needed. To specify filters to be used, simply list them as key:threshold pairs separated by comma(,). For filters that don't require threshold, leave blank after ':'.

Not all filters are available for all variant callers. Which filters can be used with your data depends on the variant caller that was used to create your VCF file.

=============  ==========================================
Filter         Remark
=============  ==========================================
Quality score  GATK and samtools mpileup
AD ratio       GATK and samtools mpileup version >=1.3
DP4 ratio      samtools mpileup version <=1.2
MQ score       GATK and samtools mpileup
MQ0 ratio      GATK only
MQ0F ratio     samtools mpileup only
GQ score       GATK and samtools mpileup
Minimum depth  GATK and samtools mpileup
Uncall GT      GATK and samtools mpileup
=============  ==========================================


If you want to filter a VCF file that was not created with either GATK or samtools, please refer to the documentation of this tool to check which data is available in your VCF files.
