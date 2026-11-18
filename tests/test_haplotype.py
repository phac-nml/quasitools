from quasitools.haplotype import Haplotype


def test_haplotype_comparison_with_other_types():
    haplotype = Haplotype('AC')

    assert haplotype == Haplotype('AC')
    assert haplotype != Haplotype('GT')
    assert not (haplotype != 'AC')
    assert not (haplotype == 'AC')
