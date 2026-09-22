"""Exact invariants of the codon-usage indices."""
import pytest

from codonamr.genetic_code import CODON2AA, DEGENERATE, FAMILY
from codonamr.metrics import (aa_composition, build_w, cai, codon_counts, enc,
                              fop, gc3s, gc_content, rscu, scuo)

from conftest import one_codon_per_aa, random_cds, uniform_usage


def test_enc_is_20_for_one_codon_per_amino_acid():
    """Wright's Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6, and every F is 1 here."""
    assert enc(one_codon_per_aa()) == pytest.approx(20.0, abs=1e-9)


def test_enc_is_near_61_for_uniform_synonymous_usage():
    """Even usage gives F = 1/k per family, so Nc approaches 61 from above.

    The finite-sample estimator slightly exceeds 61 at any achievable sequence
    length, which is why the implementation caps it, following CodonW. The test
    asserts the value is in the top of the range rather than exactly 61, so it
    would still fail if the estimator were wrong in either direction.
    """
    value = enc(uniform_usage(repeats=60))
    assert 58.0 <= value <= 61.0


def test_enc_orders_biased_below_unbiased(biased_seq, even_seq):
    assert enc(biased_seq) < enc(even_seq)


def test_enc_returns_none_when_uninformative():
    assert enc("ATGAAATTTTAA"[:9]) is None


def test_rscu_sums_to_family_size_within_each_family(biased_seq):
    table = rscu(biased_seq)
    counts = codon_counts(biased_seq)
    for aa, fam in FAMILY.items():
        if sum(counts[c] for c in fam) == 0:
            assert all(table[c] is None for c in fam)
            continue
        assert sum(table[c] for c in fam) == pytest.approx(len(fam), abs=1e-9)


def test_rscu_is_one_everywhere_under_even_usage():
    table = rscu(uniform_usage(repeats=10))
    for aa, fam in DEGENERATE.items():
        for c in fam:
            assert table[c] == pytest.approx(1.0, abs=1e-9)


def test_rscu_empty_family_is_none_not_zero():
    table = rscu("ATGATGATG")
    assert table["ATG"] == pytest.approx(1.0)
    assert table["GCT"] is None


def test_cai_of_a_maximally_biased_reference_set_against_itself_is_one():
    """w is f/f_max within each family, so a gene using only the top codon scores 1."""
    ref = [one_codon_per_aa(repeats=r) for r in (20, 30, 40)]
    w = build_w(ref)
    for seq in ref:
        assert cai(seq, w) == pytest.approx(1.0, abs=1e-12)


def test_cai_of_a_realistic_reference_set_beats_an_unrelated_gene(rng):
    """A reference set scores itself high and an evenly-using gene low."""
    ref = [random_cds(rng, 300, bias=0.9, stop="") for _ in range(30)]
    w = build_w(ref)
    ref_cai = sum(cai(s, w) for s in ref) / len(ref)
    other = cai(random_cds(rng, 300, bias=0.0, stop=""), w)
    assert ref_cai > 0.75
    assert other < 0.2
    assert ref_cai > other + 0.5


def test_cai_is_a_geometric_mean():
    w = {c: 0.5 for c in CODON2AA}
    seq = "ATGGCTGCCGCAGCGTAA"
    assert cai(seq, w) == pytest.approx(0.5)


def test_gc_content_and_gc3s():
    assert gc_content("GGCC") == 1.0
    assert gc_content("ATAT") == 0.0
    # Met and Trp third positions are excluded from GC3s by definition
    assert gc3s("ATGTGG") is None
    assert gc3s("GCGGCG") == 1.0


def test_scuo_is_zero_for_even_usage_and_one_for_single_codon_usage():
    assert scuo(uniform_usage(repeats=10)) == pytest.approx(0.0, abs=1e-9)
    assert scuo(one_codon_per_aa(repeats=10)) == pytest.approx(1.0, abs=1e-9)


def test_fop_counts_only_degenerate_codons():
    optimal = {"GCT"}
    # two Ala codons, one optimal, plus Met and Trp which must not count
    assert fop("ATGGCTGCCTGG", optimal) == pytest.approx(0.5)


def test_codon_counts_cover_all_sense_codons_and_drop_stops():
    counts = codon_counts("ATGTAA")
    assert len(counts) == 61
    assert "TAA" not in counts
    assert counts["ATG"] == 1


def test_aa_composition_matches_codon_counts(biased_seq):
    aas = aa_composition(biased_seq)
    counts = codon_counts(biased_seq)
    for aa, fam in FAMILY.items():
        assert aas.get(aa, 0) == sum(counts[c] for c in fam)
