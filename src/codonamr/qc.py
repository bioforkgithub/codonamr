"""Coding-sequence quality control.

Every metric in this package assumes a clean, in-frame CDS. Garbage in one
sequence quietly distorts a genome-level reference set, so QC is applied once,
here, and the reasons for rejection are reported rather than swallowed.
"""
from .genetic_code import STOPS, codon_list

#: why a sequence was rejected; returned by :func:`check_cds`
REASONS = ("ok", "too_short", "not_multiple_of_three", "ambiguous_bases",
           "no_stop_codon", "internal_stop", "bad_start")


def check_cds(seq, min_len=90, require_stop=True, require_start=True,
              starts=("ATG", "GTG", "TTG")):
    """Return ``(cleaned_sequence_or_None, reason)``.

    The cleaned sequence is uppercase DNA with the terminal stop removed.

    ``min_len`` defaults to 90 nt (30 codons). Note that codon-usage indices are
    unstable on short sequences: ENC in particular needs enough codons per
    synonymous family to estimate homozygosity, which is why :func:`codonamr.
    metrics.enc` returns ``None`` below its information threshold rather than a
    misleadingly precise number.
    """
    s = str(seq).upper().replace("U", "T")
    if len(s) < min_len:
        return None, "too_short"
    if len(s) % 3:
        return None, "not_multiple_of_three"
    if set(s) - set("ACGT"):
        return None, "ambiguous_bases"
    codons = codon_list(s)
    if codons[-1] in STOPS:
        codons = codons[:-1]
    elif require_stop:
        return None, "no_stop_codon"
    if any(c in STOPS for c in codons):
        return None, "internal_stop"
    if require_start and codons[0] not in starts:
        return None, "bad_start"
    return "".join(codons), "ok"


def clean_cds(seq, **kw):
    """Backwards-compatible wrapper returning only the sequence, or ``None``."""
    return check_cds(seq, **kw)[0]
