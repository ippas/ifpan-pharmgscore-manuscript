import hail as hl


class VCFFilter:
    def __init__(self):
        pass

    def min_mean_read_depth(self, mt, ratio=7):
        mt = mt.filter_rows(
            hl.agg.mean(mt.DP) >= ratio
        )
        return mt

    def min_variant_missingness(self, mt, ratio=0.1):
        if not hasattr(mt, 'variant_qc'):
            raise ValueError('variant_qc not found in mt')
        return mt.filter_rows(mt.variant_qc.call_rate >= ratio)

    def min_allele_balance(self, mt, sample_nr=1, sample_ratio=None, ratio=0.15):
        if not hasattr(mt, 'was_split'):
            raise ValueError('was_split not found in mt')
        mt = mt.filter_rows(
            hl.agg.any(
                mt.GT.is_het()
                & (hl.min(mt.AD) / hl.sum(mt.AD) >= ratio)
            )
            | hl.agg.all(~mt.GT.is_het())
        )
        return mt

    def is_indel(self):
        pass

    def hardy_weinberg(self):
        pass
