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

    name = "gatk"

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
        ref_name, _ = os.path.splitext(ref)
        opts = {"ref": ref,
                "ref_name": ref_name,
                "bam": bam,
                "picard_tools_path": os.path.join(os.environ["PICARD_TOOLS_PATH"], "CreateSequenceDictionary.jar"),
                "all_variants_file": os.path.join(variant_dir, "variants.vcf")}

        os.system("samtools faidx %(ref)s" % opts)

        os.system("java -jar %(picard_tools_path)s R=%(ref)s O=%(ref_name)s.dict" % opts)

        #  call variants
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T UnifiedGenotyper -R  %(ref)s --sample_ploidy 2 --genotype_likelihoods_model BOTH -rf BadCigar -out_mode EMIT_ALL_SITES -I %(bam)s -o %(all_variants_file)s" % opts)

        opts["non_snp_vcf"] = os.path.join(variant_dir, "variants.non_snp.vcf")
        #  select non-snps
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T SelectVariants -R %(ref)s -V %(all_variants_file)s -o %(non_snp_vcf)s  -selectType SYMBOLIC -selectType MNP -selectType INDEL -selectType MIXED -selectType NO_VARIATION" % opts)

        opts["snp_vcf"] = os.path.join(variant_dir, "variants.snp.vcf")
        # select snps
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T SelectVariants -R %(ref)s -V %(all_variants_file)s -o %(snp_vcf)s -selectType SNP" % opts)

        return opts["snp_vcf"]
