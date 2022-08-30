import hail as hl


def liftover_vcf(vcf, from_rg, to_rg, chain_file):
    rg_from = hl.get_reference(from_rg)
    rg_to = hl.get_reference(to_rg)
    if not rg_from.has_liftover(rg_to):
        rg_from.add_liftover(chain_file, rg_to)

    # https://hail.is/docs/0.2/guides/genetics.html#liftover-variants-from-one-coordinate-system-to-another
    vcf_new = vcf.annotate(
        new_locus=hl.liftover(vcf.locus, rg_to, include_strand=True),
        old_locus=vcf.locus
    )
    vcf_new = vcf_new.filter(
        hl.is_defined(vcf_new.new_locus)
        & ~vcf_new.new_locus.is_negative_strand
    )
    vcf_new = vcf_new.key_by(
        locus=vcf_new.new_locus.result,
        alleles=vcf_new.alleles
    )
    vcf_new = vcf_new.drop('new_locus')

    return vcf_new