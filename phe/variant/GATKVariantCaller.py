'''
Created on 22 Sep 2015

@author: alex
'''
import os

from phe.variant import VariantCaller


class GATKVariantCaller(VariantCaller):
    '''
    classdocs
    '''

    _default_options = None
    _threads = 1

    def __init__(self, cmd_options=None):
        '''
        Constructor
        '''
        super(GATKVariantCaller, self).__init__(cmd_options=cmd_options)


    def make_vcf(self, *args, **kwargs):

        ref = kwargs.get("ref")
        bam = kwargs.get("bam")
        variant_dir = os.path.join(kwargs.get("out_dir"), ".")

        opts = {"ref": ref,
                "bam": bam,
                "all_variants_file": os.path.join(variant_dir, "variants.vcf")}

        #  call variants
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T UnifiedGenotyper -R  %(ref)s --sample_ploidy 2 --genotype_likelihoods_model BOTH -rf BadCigar -out_mode EMIT_ALL_SITES -I %(bam)s -o %(all_variants_file)s" % opts)

        opts["non_snp_vcf"] = os.path.join(variant_dir, "variants.non_snp.vcf")
        #  select non-snps
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T SelectVariants -R %(ref)s -V %(all_variants_file)s -o %(non_snp_vcf)s  -selectType SYMBOLIC -selectType MNP -selectType INDEL -selectType MIXED -selectType NO_VARIATION" % opts)

        opts["snp_vcf"] = os.path.join(variant_dir, "variants.snp.vcf")
        # select snps
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T SelectVariants -R %(ref)s -V %(all_variants_file)s -o %(snp_vcf)s -selectType SNP" % opts)

        return True
