===================
Converting to FASTA
===================

A lot of downstream applications take on FASTA formated files, not VCFs. We have included a script for converting VCF data to FASTA format.

.. code:: bash

   $ vcfs2fasta -d <path to directory of VCFs> -o <path to output FASTA>

This tool is also able to filter out samples and columns that may be bad quality. E.g. **--sample-Ns** specifies maximum fraction of Ns present in a sample. **--Ns** specifies maximum fraction of Ns allowed per column. **--with-mixtures** can specify if mixed position should be output as over certain fraction. First, samples are sifted out then columns are checked. **--with-stats** allows to output genomic positions of the called SNPs. Each line coresponds to a character position in the FASTA file. E.g. to get the position of the 23rd snp, go to 23rd line in the positions file. With this option, if numpy and matplotlib are installed, **plots** directory will be created and a summary plot is generated, summarising called SNPs, Ns, mixtures, and bases with 0 depth. Providing **--with-dist-matrix** will also write distance matrix to the specified path.