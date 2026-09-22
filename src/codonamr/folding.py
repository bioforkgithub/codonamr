"""Local mRNA folding free energy.

ViennaRNA is an optional dependency. Import it lazily so that the codon-usage
half of the package works without it.
"""
__all__ = ["dg_windowed", "available"]

_RNA = None


def available():
    """True if ViennaRNA can be imported."""
    try:
        import RNA  # noqa: F401
        return True
    except ImportError:
        return False


def _rna():
    global _RNA
    if _RNA is None:
        try:
            import RNA
        except ImportError as exc:                       # pragma: no cover
            raise ImportError(
                "ViennaRNA is required for folding metrics. Install it with "
                "`pip install ViennaRNA`, or run codonamr with --no-folding."
            ) from exc
        _RNA = RNA
    return _RNA


def dg_windowed(seq, window=40, step=1, temperature=None):
    """Mean local folding free energy over sliding windows.

    Folds every ``window``-nt subsequence at ``step`` offsets and averages the
    minimum free energies. Window 40 / step 1 follows Victor et al. (2019).

    ``temperature`` in degrees Celsius sets the ViennaRNA model temperature. Be
    careful with thermophiles: at 80 C the model retains almost no structure, so
    the statistic loses most of its resolution and comparisons against mesophiles
    folded at 37 C are not meaningful.
    """
    RNA = _rna()
    rna = seq.upper().replace("T", "U")
    if len(rna) < window:
        return None
    if temperature is None:
        fold = RNA.fold
    else:
        md = RNA.md()
        md.temperature = float(temperature)
        def fold(s):
            fc = RNA.fold_compound(s, md)
            return fc.mfe()
    total = n = 0
    for i in range(0, len(rna) - window + 1, step):
        total += fold(rna[i:i + window])[1]
        n += 1
    return total / n if n else None
