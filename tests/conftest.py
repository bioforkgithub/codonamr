"""Fixtures: synthetic coding sequences with known codon usage.

Synthetic sequences are used because the invariants being tested are exact.
A sequence built to use one codon per amino acid has an ENC of exactly 20, and
a test that asserts that is checking the implementation rather than checking
that a downloaded sequence still looks the way it did last year.
"""
import random

import pytest

from codonamr.genetic_code import DEGENERATE, FAMILY

STOP = "TAA"


def one_codon_per_aa(repeats=40):
    """Every amino acid encoded by a single codon. ENC must be exactly 20."""
    body = "".join(FAMILY[aa][0] * repeats for aa in sorted(FAMILY)
                   if aa != "*")
    return "ATG" + body + STOP


def uniform_usage(repeats=50):
    """Every codon of every family used equally often. ENC must be near 61."""
    body = ""
    for aa in sorted(DEGENERATE):
        body += "".join(c * repeats for c in DEGENERATE[aa])
    return "ATG" + body + STOP


def random_cds(rng, n_codons=300, bias=None, start="ATG", stop=STOP):
    """A random in-frame CDS with no internal stops.

    ``bias`` in [0, 1] is the probability of taking the family's first codon,
    so 1 gives maximal bias and 0 gives even usage.
    """
    aas = [a for a in sorted(FAMILY) if a != "*"]
    out = []
    for _ in range(n_codons):
        aa = rng.choice(aas)
        fam = FAMILY[aa]
        if bias is not None and rng.random() < bias:
            out.append(fam[0])
        else:
            out.append(fam[rng.randrange(len(fam))])
    return start + "".join(out) + stop


@pytest.fixture
def rng():
    return random.Random(20240101)


@pytest.fixture
def biased_seq(rng):
    """QC-cleaned, that is no terminal stop.

    This is what the rest of the package passes around: :func:`codonamr.qc.
    check_cds` strips the stop once, at the start, and the randomisation
    designs assume it is gone because a stop codon has no synonymous family to
    permute within.
    """
    return random_cds(rng, 400, bias=0.85, stop="")


@pytest.fixture
def even_seq(rng):
    return random_cds(rng, 400, bias=0.0, stop="")


@pytest.fixture
def host_cds(rng):
    """A small synthetic genome: 300 genes with a consistent codon preference."""
    return [random_cds(rng, rng.randrange(80, 500), bias=0.55, stop="")
            for _ in range(300)]
