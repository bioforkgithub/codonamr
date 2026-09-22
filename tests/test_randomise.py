"""What each randomisation design holds fixed, asserted exactly."""
import random

import pytest

from codonamr.genetic_code import CODON2AA, FAMILY, codon_list
from codonamr.metrics import build_w, cai, codon_counts, enc, gc3s
from codonamr.randomise import (junction_counts, randomise_replace,
                                randomise_shuffle, randomise_shuffle_dinuc)

from conftest import random_cds


def _translate(seq):
    return "".join(CODON2AA[c] for c in codon_list(seq))


def _l1(a, b):
    return sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))


def test_shuffle_preserves_codon_counts_enc_cai_and_gc3s(biased_seq, rng):
    w = build_w([biased_seq])
    before = (codon_counts(biased_seq), enc(biased_seq),
              cai(biased_seq, w), gc3s(biased_seq))
    for _ in range(20):
        s = randomise_shuffle(biased_seq, rng)
        assert codon_counts(s) == before[0]
        assert enc(s) == pytest.approx(before[1], abs=1e-12)
        assert cai(s, w) == pytest.approx(before[2], abs=1e-12)
        assert gc3s(s) == pytest.approx(before[3], abs=1e-12)


def test_shuffle_preserves_the_peptide(biased_seq, rng):
    for _ in range(10):
        assert _translate(randomise_shuffle(biased_seq, rng)) == _translate(biased_seq)


def test_shuffle_actually_moves_something(biased_seq, rng):
    """A null that never changes the sequence is not a null."""
    changed = sum(randomise_shuffle(biased_seq, rng) != biased_seq
                  for _ in range(10))
    assert changed >= 9


def test_shuffle_dinuc_preserves_codon_counts_and_reports_a_residual(biased_seq, rng):
    counts = codon_counts(biased_seq)
    for _ in range(5):
        s, residual = randomise_shuffle_dinuc(biased_seq, rng)
        assert codon_counts(s) == counts
        assert _translate(s) == _translate(biased_seq)
        assert isinstance(residual, (int, float))
        assert residual >= 0
        native = junction_counts(codon_list(biased_seq))
        assert residual == _l1(native, junction_counts(codon_list(s)))


def test_shuffle_dinuc_is_never_worse_than_a_plain_shuffle(biased_seq):
    """The dinucleotide design descends from a plain shuffle, so it cannot lose."""
    native = junction_counts(codon_list(biased_seq))
    plain = sum(_l1(native, junction_counts(codon_list(
        randomise_shuffle(biased_seq, random.Random(i)))))
        for i in range(10)) / 10.0
    controlled = sum(randomise_shuffle_dinuc(biased_seq, random.Random(i))[1]
                     for i in range(10)) / 10.0
    assert controlled <= plain


def test_replace_preserves_the_peptide_and_changes_composition(biased_seq, rng):
    preferred = {aa: fam[0] for aa, fam in FAMILY.items()}
    maximal = randomise_replace(biased_seq, preferred, 1.0, rng)
    assert _translate(maximal) == _translate(biased_seq)
    assert enc(maximal) == pytest.approx(20.0, abs=1e-9)


def test_replace_at_p_equals_one_over_k_is_roughly_uniform(rng):
    seq = random_cds(rng, 4000, bias=0.95, stop="")
    preferred = {aa: fam[0] for aa, fam in FAMILY.items()}
    uniform = randomise_replace(seq, preferred, 0.0, rng)
    assert enc(uniform) > enc(seq) + 10


def test_junction_counts_only_spans_codon_boundaries():
    j = junction_counts(["ATG", "GCT", "TAA"])
    assert j == {"GG": 1, "TT": 1}
    assert sum(j.values()) == 2
