"""The length behaviour that :func:`codonamr.host.cufs` warns about.

The warning in ``host.cufs`` is a quantitative claim about how these indices
behave when gene length varies, and a published analysis was nearly reversed
by getting it wrong. Two things are pinned here.

1. The behaviour itself. With codon preference held exactly fixed, shortening
   a gene inflates its CUFS and RSCU distance to the reference it was drawn
   from, while leaving CAI unbiased. That asymmetry is the whole content of
   the warning.
2. The sign printed in the docstring. An earlier version of the warning gave
   CAI's length correlation as -0.20 when the measured value is +0.20, which
   tells a reader the bias runs the same way as the distances when it runs the
   other way. Regression test, since a docstring is the only place most users
   will meet this.
"""
import random

import pytest

from codonamr.genetic_code import CODON2AA, FAMILY
from codonamr.host import cufs, rscu_distance
from codonamr.metrics import build_w, cai, rscu

DEGENERATE_CODONS = [c for c, a in CODON2AA.items()
                     if a != "*" and len(FAMILY[a]) > 1]


def _spearman(x, y):
    """Spearman's rho without pulling scipy into a unit test."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else 0.0


def _corpus(seed=7, n_codons=60000):
    """A long sequence with a fixed, biased synonymous codon preference."""
    rnd = random.Random(seed)
    weights = {}
    for aa, fam in FAMILY.items():
        ws = [rnd.random() ** 2 + 0.05 for _ in fam]
        tot = sum(ws)
        for c, w in zip(fam, ws):
            weights[c] = w / tot
    codons, aas = [], [a for a in FAMILY if a != "*"]
    for _ in range(n_codons):
        aa = rnd.choice(aas)
        fam = FAMILY[aa]
        codons.append(rnd.choices(fam, [weights[c] for c in fam])[0])
    return codons


LENGTHS = [60, 90, 130, 190, 280, 400, 600, 900, 1400, 2000]


@pytest.fixture(scope="module")
def sampled():
    """Sub-genes of many lengths drawn from ONE codon distribution.

    Every sub-gene has the same expected codon preference, so any systematic
    relationship between length and a score is an artefact of the score.
    """
    corpus = _corpus()
    reference = rscu(dict((c, corpus.count(c)) for c in set(corpus)))
    w = build_w(["".join(corpus)])
    rnd = random.Random(11)
    out = []
    for n in LENGTHS:
        for _ in range(25):
            i = rnd.randrange(0, len(corpus) - n)
            seq = "".join(corpus[i:i + n])
            out.append((n, seq, rscu(seq)))
    return reference, w, out


def test_cufs_is_inflated_by_short_length(sampled):
    reference, _w, rows = sampled
    lens = [n for n, _s, _r in rows]
    vals = [cufs(r, reference) for _n, _s, r in rows]
    rho = _spearman(lens, vals)
    assert rho < -0.5, (
        "CUFS to the distribution a gene was drawn from should fall with "
        "length even though preference is constant; measured rho %.3f" % rho)


def test_rscu_distance_is_inflated_by_short_length(sampled):
    reference, _w, rows = sampled
    lens = [n for n, _s, _r in rows]
    for metric, bound in (("euclidean", -0.5), ("correlation", -0.5)):
        vals = [rscu_distance(r, reference, metric=metric)
                for _n, _s, r in rows]
        rho = _spearman(lens, vals)
        assert rho < bound, "%s rho %.3f" % (metric, rho)


def test_cai_is_not_inflated_by_short_length(sampled):
    """CAI is a per-codon mean, so it carries no systematic length bias.

    This is what makes it usable across genes of different length and it is
    the half of the warning that the wrong sign contradicted.
    """
    _reference, w, rows = sampled
    lens = [n for n, _s, _r in rows]
    vals = [cai(s, w) for _n, s, _r in rows]
    rho = _spearman(lens, vals)
    assert abs(rho) < 0.25, (
        "CAI should show no systematic drift with length at fixed codon "
        "preference; measured rho %.3f" % rho)


def test_short_genes_are_compared_on_fewer_families(sampled):
    """The countable half of the mechanism, quoted in the cufs warning."""
    _reference, _w, rows = sampled
    lens = [n for n, _s, _r in rows]
    n_fam = []
    for _n, _s, r in rows:
        n_fam.append(sum(1 for aa, fam in FAMILY.items()
                         if len(fam) > 1
                         and all(r.get(c) is not None for c in fam)))
    assert _spearman(lens, n_fam) > 0.4
    assert min(n_fam) < max(n_fam)


def test_cufs_warning_gives_cai_a_positive_length_correlation():
    """Regression: the published warning once printed this sign backwards."""
    doc = cufs.__doc__
    assert "for CAI" in doc
    i = doc.index("for CAI")
    window = doc[max(0, i - 30):i]
    assert "+0.20" in window, (
        "codonamr.host.cufs must quote CAI's length correlation as +0.20. It "
        "is positive: long genes score slightly higher CAI, the opposite of "
        "what happens to CUFS and RSCU distance. Window was %r" % window)
    assert "-0.20 for CAI" not in doc
