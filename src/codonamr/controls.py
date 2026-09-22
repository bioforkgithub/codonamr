"""Controls, and the validation of those controls.

Most reported associations between codon usage and anything else are
confounded. Two confounders dominate. The first is length: ENC and every other
family-based index is a small-sample statistic, so a short gene has a biased and
noisy value, and resistance genes are short. Comparing a 150-codon gene against
a genome mean dominated by 300-codon genes produces a difference that is an
artefact of length alone. The second is composition: GC content drives codon
choice mutationally, and controlling for GC as one number leaves the strand-
asymmetric part of composition uncontrolled, so an association with, say, A
content survives as an apparent association with codon usage.

The functions here build those controls, and then check that they work. A
control nobody has tested is an assumption. :func:`injection_recovery` injects
an association of known size and asks what comes back; :func:`false_positive_rate`
injects none and asks how often something comes back anyway. Both are meant to
be run and reported alongside the analysis, which is what ``codonamr validate``
does.

References
----------
ENC small-sample behaviour  Fuglsang A (2004) Biochem Biophys Res Commun 317:957-964.
Partial correlation         Kendall MG, Stuart A (1973) The Advanced Theory of Statistics, vol 2.
Permutation and empirical p North BV, Curtis D, Sham PC (2002) Am J Hum Genet 71:439-441.
"""
import math
import random
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from .genetic_code import codon_list

__all__ = ["length_matched_null", "LengthMatchedNull", "empirical_z",
           "base_fractions", "spearman_partial", "composition_control",
           "injection_recovery", "false_positive_rate"]


# ------------------------------------------------------- length matching
@dataclass
class LengthMatchedNull:
    """Host genes matched to a query gene's length, with the matching recorded.

    Iterating over this yields the sequences, so it can be used directly where
    a list of sequences is expected, while the diagnostics stay attached.
    """
    sequences: list = field(default_factory=list)
    target_codons: int = 0
    pool_size: int = 0
    tolerance: float = 0.0
    with_replacement: bool = False
    note: str = ""

    def __iter__(self):
        return iter(self.sequences)

    def __len__(self):
        return len(self.sequences)


def length_matched_null(seq, host_cds, n=100, rng=None, tolerance=0.2,
                        min_pool=20, widen=(0.2, 0.35, 0.5, 1.0)):
    """Draw host genes of the same length as ``seq``.

    Parameters
    ----------
    seq : str
        The query coding sequence.
    host_cds : sequence of str
        The host's coding sequences, already QC-passed.
    n : int
        How many null genes to draw.
    rng : random.Random, optional
        Seeded generator. One is created from the system entropy if omitted,
        which makes the result irreproducible, so the pipeline always passes one.
    tolerance : float
        Starting relative tolerance on codon length, so 0.2 accepts genes
        between 0.8 and 1.2 times the query's length.
    min_pool : int
        Fewest genes the matched pool may contain before the tolerance is
        widened.
    widen : tuple of float
        Tolerances tried in order when the pool is too small.

    Returns
    -------
    LengthMatchedNull

    Notes
    -----
    ENC is estimated from the homozygosity of each synonymous family, and that
    estimate is biased upward and highly variable when a family holds only a few
    codons. Below roughly 200 codons the bias is large enough to dominate any
    biological signal, and below about 100 codons ENC is close to uninformative.
    The consequence is that a short gene will look less biased than a long one
    from the same genome, with no difference in selection at all. Comparing the
    query against host genes of its own length removes that, because the bias
    applies equally to both sides. It is the reason this function exists, and it
    is why the pipeline reports a length-matched z score rather than a raw
    difference from the genome mean.

    If the pool is smaller than ``n`` after widening, genes are drawn with
    replacement and the fact is recorded in ``with_replacement``. A null built
    by resampling twelve genes a hundred times is not a hundred independent
    draws, and any p-value from it should be read as a lower bound on the true
    one.
    """
    rng = rng or random.Random()
    target = len(codon_list(seq))
    pool, used_tol, note = [], tolerance, ""
    for tol in (t for t in (tolerance,) + tuple(widen) if t >= tolerance):
        lo, hi = target * (1 - tol), target * (1 + tol)
        pool = [s for s in host_cds if lo <= len(codon_list(s)) <= hi]
        used_tol = tol
        if len(pool) >= max(min_pool, 1):
            break
    if len(pool) < max(min_pool, 1):
        ranked = sorted(host_cds,
                        key=lambda s: abs(math.log((len(codon_list(s)) or 1)
                                                   / max(target, 1))))
        pool = ranked[:max(min_pool, n)]
        used_tol = float("inf")
        note = ("no gene pool within the widest tolerance, fell back to the "
                "nearest genes by length ratio")
    if not pool:
        return LengthMatchedNull([], target, 0, used_tol, False,
                                 "no host sequences available")
    replace = len(pool) < n
    if replace:
        draw = [pool[rng.randrange(len(pool))] for _ in range(n)]
    else:
        draw = rng.sample(pool, n)
    return LengthMatchedNull(draw, target, len(pool), used_tol, replace, note)


def empirical_z(value, null_values):
    """Standardise an observation against a null sample, with an empirical p.

    The p-value uses the ``(r + 1) / (m + 1)`` correction of North, Curtis and
    Sham (2002), so it can never be reported as exactly zero from a finite
    number of draws. With 100 permutations the smallest attainable two-sided p
    is about 0.02, which is a property of the design and not of the data.
    """
    vals = np.asarray([v for v in null_values if v is not None], dtype=float)
    vals = vals[np.isfinite(vals)]
    if value is None or not np.isfinite(value) or vals.size < 2:
        return {"z": None, "p": None, "n_null": int(vals.size),
                "null_mean": None, "null_sd": None}
    mu, sd = float(vals.mean()), float(vals.std(ddof=1))
    z = (value - mu) / sd if sd > 0 else None
    more = int(np.sum(np.abs(vals - mu) >= abs(value - mu)))
    p = (more + 1) / (vals.size + 1)
    return {"z": z, "p": min(p, 1.0), "n_null": int(vals.size),
            "null_mean": mu, "null_sd": sd}


# ------------------------------------------------------ composition control
def base_fractions(seq):
    """Fractions of A, C, G and T, as a 4-tuple, over the whole sequence."""
    s = str(seq).upper()
    n = len(s)
    if not n:
        return (None, None, None, None)
    return tuple(s.count(b) / n for b in "ACGT")


def _clean(x, y, covariates):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    c = np.asarray(covariates, dtype=float)
    if c.ndim == 1:
        c = c[:, None]
    ok = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(c), axis=1)
    return x[ok], y[ok], c[ok]


def _residuals(v, design):
    beta, *_ = np.linalg.lstsq(design, v, rcond=None)
    return v - design @ beta


def spearman_partial(x, y, covariates):
    """Spearman partial correlation of ``x`` and ``y`` given ``covariates``.

    Everything is rank-transformed first, then ``x`` and ``y`` are regressed on
    the ranked covariates and Pearson's correlation is taken between the
    residuals. This is the usual rank analogue of a partial correlation. It
    removes monotone dependence on the covariates, not arbitrary dependence: a
    U-shaped relationship with GC survives the control and can still generate a
    spurious partial correlation.

    Parameters
    ----------
    x, y : array-like
    covariates : array-like
        One column per covariate, or a 1-D array for a single covariate. Rows
        with a non-finite value anywhere are dropped pairwise-complete.

    Returns
    -------
    dict
        ``r``, ``p``, ``n``, ``k`` (number of covariates) and ``df``.
        ``r`` is ``None`` when fewer than ``k + 4`` complete rows remain or a
        residual has no variance.

    Notes
    -----
    The p-value comes from ``t = r * sqrt(df / (1 - r^2))`` with
    ``df = n - 2 - k``. Because the ranks are used as if they were the original
    values, this p is approximate; it is well behaved for n above about 30 and
    should not be trusted for very small samples. Collinear covariates are
    handled by a least-squares solve with rank detection rather than an inverse,
    so a redundant column reduces the effective degrees of freedom instead of
    raising an error.
    """
    x, y, c = _clean(x, y, covariates)
    n, k = x.size, c.shape[1]
    if n < k + 4:
        return {"r": None, "p": None, "n": int(n), "k": int(k), "df": 0}
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rc = np.column_stack([stats.rankdata(c[:, j]) for j in range(k)])
    design = np.column_stack([np.ones(n), rc])
    rank = int(np.linalg.matrix_rank(design))
    ex, ey = _residuals(rx, design), _residuals(ry, design)
    if ex.std() == 0 or ey.std() == 0:
        return {"r": None, "p": None, "n": int(n), "k": int(k), "df": 0}
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - rank - 1
    if df <= 0:
        return {"r": r, "p": None, "n": int(n), "k": int(k), "df": 0}
    r_clipped = min(max(r, -0.999999999), 0.999999999)
    t = r_clipped * math.sqrt(df / (1 - r_clipped ** 2))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"r": r, "p": p, "n": int(n), "k": int(k), "df": int(df)}


def composition_control(values, gc_per_nucleotide, response, method="spearman",
                        drop_base=3):
    """Correlate ``values`` with ``response``, controlling for base composition.

    Parameters
    ----------
    values : array-like
        The codon-usage quantity under test, one value per gene, for example
        ENC or CAI.
    gc_per_nucleotide : array-like, shape (n, 4)
        Per-nucleotide composition, columns in A, C, G, T order, as returned by
        :func:`base_fractions`. A 1-D array is accepted and treated as a single
        composition covariate, which reproduces the usual GC-only control for
        comparison.
    response : array-like
        The outcome being related to codon usage, one value per gene.
    method : {"spearman", "pearson"}
        ``"pearson"`` skips the rank transform. Use it only when both variables
        are plausibly linear and unbounded; ENC and CAI are neither.
    drop_base : int or None
        Index of the composition column to drop, default 3 (T). The four
        fractions sum to 1, so using all four makes the design matrix singular.
        Dropping one loses nothing: the remaining three span the same space.
        Set to ``None`` to keep all four and rely on the rank-deficient solve.

    Returns
    -------
    dict
        ``r_raw`` and ``p_raw`` for the uncontrolled association, ``r_partial``
        and ``p_partial`` for the controlled one, plus ``n``, ``k`` and ``df``.

    Notes
    -----
    Controlling for GC alone treats A and T, and C and G, as interchangeable.
    They are not: replication-associated strand asymmetry, amino-acid
    composition and transcription-coupled repair all act on individual bases,
    and a mobile gene acquired from a compositionally different genome
    typically differs in more than its GC total. Regressing on A, C and G
    fractions separately therefore removes strictly more confounding than a GC
    control, at the cost of three degrees of freedom. If the partial
    correlation survives that and the GC-only one does not, it is the GC-only
    result that was wrong.
    """
    comp = np.asarray(gc_per_nucleotide, dtype=float)
    if comp.ndim == 1:
        comp = comp[:, None]
    if drop_base is not None and comp.shape[1] == 4:
        comp = np.delete(comp, drop_base, axis=1)
    x = np.asarray(values, dtype=float)
    y = np.asarray(response, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(comp), axis=1)
    xs, ys = x[ok], y[ok]
    if xs.size < 4:
        return {"r_raw": None, "p_raw": None, "r_partial": None,
                "p_partial": None, "n": int(xs.size), "k": int(comp.shape[1]),
                "df": 0}
    if method == "spearman":
        raw = stats.spearmanr(xs, ys)
        part = spearman_partial(xs, ys, comp[ok])
    elif method == "pearson":
        raw = stats.pearsonr(xs, ys)
        part = spearman_partial(xs, ys, comp[ok])
        part["note"] = "partial term always uses ranks"
    else:
        raise ValueError("method must be 'spearman' or 'pearson'")
    return {"r_raw": float(raw[0]), "p_raw": float(raw[1]),
            "r_partial": part["r"], "p_partial": part["p"],
            "n": part["n"], "k": part["k"], "df": part["df"]}


# ---------------------------------------------------- validating the control
def _simulate(n, effect, confound, rng, sequences=None, skew=0.5):
    """One replicate: a confounded pair with a known partial correlation.

    Composition is either resampled from real sequences or drawn. A
    standardised confounder ``c`` is built from it as a mixture of GC content
    and AT skew,

    ``c = (1 - skew) * z(G+C) + skew * z(A-T)``

    and then

    ``x = confound * c + u`` and ``y = confound * c + rho*u + sqrt(1-rho^2)*v``

    with ``u`` and ``v`` independent standard normals. By construction the
    population partial correlation of x and y given c is exactly ``rho``, while
    the raw correlation is inflated by the shared ``confound * c`` term.

    The ``skew`` term is the point of the design. A GC-only control cannot see
    it, because AT skew leaves GC unchanged, so the GC-only row stays
    confounded in proportion to ``skew`` while the per-nucleotide row does not.
    Set ``skew=0`` to make the confounder pure GC, in which case the two
    controls should and do perform identically.
    """
    if sequences is not None:
        comp = np.asarray(sequences, dtype=float)
        if comp.ndim != 2 or comp.shape[1] != 4:
            raise ValueError("sequences must be an (m, 4) array of A,C,G,T fractions")
        comp = comp[rng.integers(0, comp.shape[0], size=n)]
    else:
        gc = rng.uniform(0.25, 0.75, size=n)
        at_skew = rng.normal(0, 0.05, size=n)
        gc_skew = rng.normal(0, 0.02, size=n)
        comp = np.column_stack([(1 - gc) / 2 + at_skew, gc / 2 + gc_skew,
                                gc / 2 - gc_skew, (1 - gc) / 2 - at_skew])
        comp = np.clip(comp, 1e-6, None)
        comp = comp / comp.sum(axis=1, keepdims=True)

    def _z(v):
        v = np.asarray(v, dtype=float)
        sd = v.std()
        return (v - v.mean()) / (sd if sd else 1.0)

    c = (1 - skew) * _z(comp[:, 1] + comp[:, 2]) + skew * _z(comp[:, 0] - comp[:, 3])
    c = _z(c)
    u = rng.normal(size=n)
    v = rng.normal(size=n)
    rho = float(effect)
    x = confound * c + u
    y = confound * c + rho * u + math.sqrt(max(1 - rho ** 2, 0.0)) * v
    return x, y, comp


def _run_replicates(n, effect, confound, n_replicates, seed, alpha, sequences,
                    skew):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(int(n_replicates)):
        x, y, comp = _simulate(n, effect, confound, rng, sequences, skew)
        gc_only = (comp[:, 1] + comp[:, 2])[:, None]
        res = composition_control(x, comp, y)
        res_gc = composition_control(x, gc_only, y)
        rows.append({
            "r_none": res["r_raw"], "p_none": res["p_raw"],
            "r_gc": res_gc["r_partial"], "p_gc": res_gc["p_partial"],
            "r_composition": res["r_partial"], "p_composition": res["p_partial"],
        })
    df = pd.DataFrame(rows)
    out = []
    labels = {"none": "no control", "gc": "GC only",
              "composition": "per-nucleotide composition"}
    for key, label in labels.items():
        r = df["r_" + key].astype(float)
        p = df["p_" + key].astype(float)
        out.append({
            "control": label,
            "n_genes": int(n),
            "replicates": int(n_replicates),
            "true_partial_r": float(effect),
            "mean_r": float(r.mean()),
            "sd_r": float(r.std(ddof=1)) if len(r) > 1 else float("nan"),
            "bias": float(r.mean() - effect),
            "rejection_rate": float((p < alpha).mean()),
        })
    return pd.DataFrame(out)


def injection_recovery(n=200, effect=0.3, confound=0.8, n_replicates=200,
                       seed=0, alpha=0.05, sequences=None, skew=0.5):
    """Inject an association of known size and report what each control returns.

    A control is only useful if it leaves real signal intact. This plants a
    partial correlation of exactly ``effect`` between two variables that also
    share a composition-driven confounder of strength ``confound``, then runs
    three analyses on the same data: no control, a GC-only control, and the
    per-nucleotide composition control.

    Parameters
    ----------
    n : int
        Genes per replicate.
    effect : float
        The true partial correlation to plant. ``0`` gives the null case, which
        is what :func:`false_positive_rate` uses.
    confound : float
        Strength of the shared dependence on composition. At 0 there is nothing
        to control for and all three rows should agree.
    n_replicates : int
        Independent replicates.
    seed : int
        Seeds a numpy Generator, so the whole table is reproducible.
    alpha : float
        Significance threshold for ``rejection_rate``, which is statistical
        power here and the false positive rate when ``effect`` is 0.
    sequences : array-like, shape (m, 4), optional
        Real per-nucleotide composition rows to resample, from
        :func:`base_fractions`. Using the composition of the genes actually
        being analysed makes the validation specific to that dataset, which is
        what the pipeline does.
    skew : float
        Fraction of the confounder carried by AT skew rather than by GC
        content. At 0 the confounder is pure GC and a GC-only control is
        sufficient; at 0.5, the default, half of it is invisible to GC.

    Returns
    -------
    pandas.DataFrame
        One row per control, with ``mean_r``, ``bias`` and ``rejection_rate``.

    Notes
    -----
    Expect the uncontrolled row to be biased upward by roughly
    ``confound**2 / (1 + confound**2)`` and the composition row to sit close to
    ``effect``. A small negative bias in the controlled row is expected and is
    not a failure: the estimate is a Spearman correlation of a Pearson-normal
    quantity, and the rank transform attenuates it by a known factor of about
    ``(6/pi) * arcsin(r/2) / r``. Judge the control by the bias relative to the
    uncontrolled row, not by an exact match to ``effect``.
    """
    return _run_replicates(n, effect, confound, n_replicates, seed, alpha,
                           sequences, skew)


def false_positive_rate(n=200, confound=0.8, n_replicates=500, seed=0,
                        alpha=0.05, sequences=None, skew=0.5):
    """Inject no association at all and report how often one is found anyway.

    Same design as :func:`injection_recovery` with ``effect=0``. The
    ``rejection_rate`` column is then the false positive rate and should be
    close to ``alpha`` for a control that works. The uncontrolled row will be
    far above ``alpha`` whenever ``confound`` is non-zero, and the size of that
    gap is the argument for controlling at all.

    With the default ``skew``, half the confounder is strand-asymmetric, so the
    GC-only row should also sit above ``alpha``. That gap between the GC-only
    and per-nucleotide rows is the concrete argument for controlling on the
    four base fractions separately rather than on GC.

    Reporting this number next to a result is the difference between a claim
    that a control was applied and evidence that it was needed.
    """
    return _run_replicates(n, 0.0, confound, n_replicates, seed, alpha,
                           sequences, skew)
