"""Figures and the Markdown report."""
import numpy as np
import pandas as pd
import pytest

from codonamr.report import (PALETTE, plot_enc_gc3s, plot_host_compatibility,
                             plot_neutrality, plot_pr2, plot_rscu_heatmap,
                             wright_expected_enc)


def test_wright_curve_hits_its_known_points():
    """At GC3s = 0.5 the expected Nc is 2 + 0.5 + 29/0.5 = 60.5."""
    assert wright_expected_enc(0.5) == pytest.approx(60.5)
    # the curve is symmetric about 0.5 apart from the linear s term
    assert wright_expected_enc(0.2) - 0.2 == pytest.approx(
        wright_expected_enc(0.8) - 0.8)
    assert wright_expected_enc(np.array([0.3, 0.7])).shape == (2,)


def test_wright_curve_never_exceeds_61_in_the_usable_range():
    s = np.linspace(0.05, 0.95, 200)
    assert wright_expected_enc(s).max() <= 61.0 + 1e-9


def _frame(n=12):
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "id": ["g%d" % i for i in range(n)],
        "gc3s": rng.uniform(0.3, 0.7, n),
        "enc": rng.uniform(35, 58, n),
        "gc12": rng.uniform(0.4, 0.6, n),
        "gc3": rng.uniform(0.3, 0.7, n),
        "pr2_a3_at3": rng.uniform(0.3, 0.7, n),
        "pr2_g3_gc3": rng.uniform(0.3, 0.7, n),
        "cai_reference": rng.uniform(0.3, 0.8, n),
        "cufs_host": rng.uniform(0.0, 0.4, n),
    })


def test_every_figure_writes_a_png_and_a_pdf(tmp_path):
    df = _frame()
    host = _frame(40)
    for fn, args in ((plot_enc_gc3s, (df, host)), (plot_neutrality, (df,)),
                     (plot_pr2, (df,)), (plot_host_compatibility, (df, host))):
        png, pdf = fn(*args, out_dir=tmp_path)
        assert png.exists() and png.stat().st_size > 5000
        assert pdf.exists() and pdf.suffix == ".pdf"


def test_rscu_heatmap_handles_missing_families(tmp_path):
    from codonamr.genetic_code import CODON2AA, CODONS
    codons = [c for c in CODONS if CODON2AA[c] != "*"]
    rng = np.random.default_rng(1)
    data = pd.DataFrame(rng.uniform(0, 3, size=(5, len(codons))),
                        columns=codons)
    data.insert(0, "id", ["g%d" % i for i in range(5)])
    data.loc[0, "GCT"] = np.nan
    png, pdf = plot_rscu_heatmap(data, {c: 1.0 for c in codons},
                                 out_dir=tmp_path)
    assert png.exists()


def test_palette_is_the_published_colourblind_safe_one():
    """Okabe and Ito 2008. Changing these values needs a fresh CVD check."""
    assert PALETTE["blue"] == "#0072B2"
    assert PALETTE["vermillion"] == "#D55E00"
    assert PALETTE["green"] == "#009E73"
    assert len({v.lower() for v in PALETTE.values()}) == len(PALETTE)


def test_plot_functions_tolerate_empty_input(tmp_path):
    empty = pd.DataFrame(columns=["id", "gc3s", "enc", "gc12", "gc3",
                                  "pr2_a3_at3", "pr2_g3_gc3",
                                  "cai_reference", "cufs_host"], dtype=float)
    png, _ = plot_enc_gc3s(empty, None, out_dir=tmp_path, name="empty")
    assert png.exists()
