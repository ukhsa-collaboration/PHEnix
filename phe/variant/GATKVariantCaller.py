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
    _default_options = "-nt 1"

    def __init__(self, cmd_options=None):
        '''
        Constructor
        '''
        if cmd_options is None:
            cmd_options = self._default_options

        super(GATKVariantCaller, self).__init__(cmd_options=cmd_options)

    def make_vcf(self, *args, **kwargs):

        ref = kwargs.get("ref")
        bam = kwargs.get("bam")

        if kwargs.get("out_dir") is None:
            kwargs["out_dir"] = "."

        variant_dir = os.path.abspath(kwargs.get("out_dir"))
        ref_name, _ = os.path.splitext(ref)
        opts = {"ref": ref,
                "ref_name": ref_name,
                "bam": bam,
                "picard_tools_path": os.path.join(os.environ["PICARD_TOOLS_PATH"], "CreateSequenceDictionary.jar"),
                "ploidy": 2,
                "glm": "BOTH",
                "all_variants_file": os.path.join(variant_dir, "variants.vcf"),
                "extra_cmd_options": self.cmd_options}

        os.system("samtools faidx %(ref)s" % opts)

        os.system("java -jar %(picard_tools_path)s R=%(ref)s O=%(ref_name)s.dict" % opts)

        #  call variants
        # FIXME: Sample ploidy = 2?
        os.system("java -XX:+UseSerialGC -jar $GATK_JAR -T UnifiedGenotyper -R  %(ref)s --sample_ploidy %(ploidy)s --genotype_likelihoods_model %(glm)s -rf BadCigar -out_mode EMIT_ALL_SITES -I %(bam)s -o %(all_variants_file)s %(extra_cmd_options)s" % opts)

        return opts["all_variants_file"]
