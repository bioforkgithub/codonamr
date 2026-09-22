"""The controls, and the validation of the controls."""
import math
import random

import numpy as np
import pytest

from codonamr.controls import (base_fractions, composition_control,
                               empirical_z, false_positive_rate,
                               injection_recovery, length_matched_null,
                               spearman_partial)
from codonamr.genetic_code import codon_list

from conftest import random_cds


def _known_partial(n, rho, confound, seed):
    """x and y share a confounder c and have partial correlation exactly rho."""
    rng = np.random.default_rng(seed)
    c = rng.normal(size=n)
    u = rng.normal(size=n)
    v = rng.normal(size=n)
    x = confound * c + u
    y = confound * c + rho * u + math.sqrt(1 - rho ** 2) * v
    return x, y, c


def test_spearman_partial_recovers_a_known_partial_correlation():
    """Planted rho = 0.5, and Spearman attenuates it by (6/pi) arcsin(rho/2)."""
    x, y, c = _known_partial(6000, 0.5, 0.9, seed=7)
    expected = (6 / math.pi) * math.asin(0.5 / 2)
    res = spearman_partial(x, y, c)
    assert res["r"] == pytest.approx(expected, abs=0.04)
    assert res["p"] < 1e-20
    assert res["n"] == 6000
    assert res["k"] == 1
    assert res["df"] == 6000 - 3


def test_spearman_partial_removes_a_pure_confounder():
    """With rho = 0, the raw correlation is large and the partial one is not."""
    x, y, c = _known_partial(4000, 0.0, 0.9, seed=11)
    from scipy import stats
    raw = stats.spearmanr(x, y)[0]
    res = spearman_partial(x, y, c)
    assert raw > 0.35
    assert abs(res["r"]) < 0.05
    assert res["p"] > 0.01


def test_spearman_partial_reports_none_when_there_is_no_data():
    res = spearman_partial([1.0, 2.0], [1.0, 2.0], [[1.0], [2.0]])
    assert res["r"] is None


def test_spearman_partial_drops_non_finite_rows():
    x, y, c = _known_partial(500, 0.4, 0.5, seed=3)
    x = x.copy()
    x[:10] = np.nan
    res = spearman_partial(x, y, c)
    assert res["n"] == 490


def test_composition_control_uses_all_four_base_fractions():
    rng = np.random.default_rng(5)
    n = 800
    comp = rng.dirichlet([8, 8, 8, 8], size=n)
    c = (comp[:, 0] - comp[:, 3])
    c = (c - c.mean()) / c.std()
    x = c + rng.normal(size=n)
    y = c + rng.normal(size=n)
    res = composition_control(x, comp, y)
    assert res["r_raw"] > 0.25
    assert abs(res["r_partial"]) < 0.12
    assert res["k"] == 3         # T dropped, three fractions span the same space


def test_base_fractions_sum_to_one():
    a, c, g, t = base_fractions("ACGTACGT")
    assert a + c + g + t == pytest.approx(1.0)
    assert a == pytest.approx(0.25)
    assert base_fractions("") == (None, None, None, None)


def test_length_matched_null_returns_genes_of_matching_length(host_cds, rng):
    query = random_cds(rng, 120)
    null = length_matched_null(query, host_cds, n=50, rng=rng, tolerance=0.2)
    target = len(codon_list(query))
    assert len(null) == 50
    assert null.target_codons == target
    for s in null:
        assert 0.8 * target <= len(codon_list(s)) <= 1.2 * target


def test_length_matched_null_is_reproducible(host_cds):
    a = length_matched_null("ATG" + "GCT" * 150, host_cds, n=20,
                            rng=random.Random(1))
    b = length_matched_null("ATG" + "GCT" * 150, host_cds, n=20,
                            rng=random.Random(1))
    assert a.sequences == b.sequences


def test_length_matched_null_records_when_it_had_to_resample(rng):
    pool = ["ATG" + "GCT" * 100] * 3
    null = length_matched_null("ATG" + "GCT" * 100, pool, n=30, rng=rng,
                               min_pool=1)
    assert null.with_replacement is True
    assert len(null) == 30


def test_length_matched_null_widens_and_says_so(rng):
    pool = [random_cds(rng, 500) for _ in range(30)]
    null = length_matched_null(random_cds(rng, 50), pool, n=10, rng=rng)
    assert null.note or null.tolerance > 0.2


def test_empirical_z_never_reports_p_equals_zero():
    res = empirical_z(100.0, list(range(20)))
    assert res["p"] == pytest.approx(1 / 21)
    assert res["z"] > 3


def test_empirical_z_handles_missing_values():
    assert empirical_z(None, [1, 2, 3])["z"] is None
    assert empirical_z(1.0, [None])["z"] is None


def test_injection_recovery_returns_the_planted_effect_only_with_the_control():
    table = injection_recovery(n=200, effect=0.3, confound=0.8, skew=0.5,
                               n_replicates=60, seed=4).set_index("control")
    assert table.loc["no control", "mean_r"] > 0.45
    assert table.loc["per-nucleotide composition", "mean_r"] == pytest.approx(
        0.3, abs=0.06)
    assert abs(table.loc["per-nucleotide composition", "bias"]) < \
        abs(table.loc["GC only", "bias"])


def test_false_positive_rate_is_near_alpha_only_with_the_control():
    table = false_positive_rate(n=200, confound=0.8, skew=0.5,
                                n_replicates=150, seed=6,
                                alpha=0.05).set_index("control")
    assert table.loc["no control", "rejection_rate"] > 0.5
    assert table.loc["GC only", "rejection_rate"] > 0.2
    assert table.loc["per-nucleotide composition", "rejection_rate"] < 0.15
    assert (table["true_partial_r"] == 0).all()


def test_validation_is_deterministic_under_a_seed():
    a = injection_recovery(n=100, n_replicates=30, seed=99)
    b = injection_recovery(n=100, n_replicates=30, seed=99)
    assert a.equals(b)


def test_gc_only_control_is_enough_when_the_confounder_is_pure_gc():
    """The per-nucleotide argument is about strand asymmetry, not about GC."""
    table = false_positive_rate(n=200, confound=0.8, skew=0.0,
                                n_replicates=150, seed=8).set_index("control")
    assert table.loc["GC only", "rejection_rate"] < 0.15
