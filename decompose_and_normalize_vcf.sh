# setup
VCF=vcf_files/fam1_variants.vcf.gz
NORMVCF=norm_vcf_files/fam1_variants_norm.vcf.gz
REF="/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"

# decompose, normalize and annotate VCF with snpEff.
# NOTE: can also swap snpEff with VEP
zless $VCF \
    | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
    | vt decompose -s - \
    | vt normalize -r $REF - \
    | snpEff GRCh37.75 \
    | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF

# # load the pre-processed VCF into GEMINI
# gemini load --cores 3 -t snpEff -v $NORMVCF $db
#
# # query away
# gemini query -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
#              --gt-filter "gt_types.mom == HET and \
#                           gt_types.dad == HET and \
#                           gt_types.kid == HOM_ALT" \
#              $db