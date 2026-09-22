"""Every QC failure mode, checked by its reason string."""
from codonamr.qc import REASONS, check_cds, clean_cds


def test_clean_sequence_passes_and_loses_its_stop():
    seq = "ATG" + "GCT" * 40 + "TAA"
    cleaned, reason = check_cds(seq)
    assert reason == "ok"
    assert cleaned == "ATG" + "GCT" * 40
    assert "TAA" not in cleaned[-3:]


def test_too_short():
    assert check_cds("ATGGCTTAA") == (None, "too_short")


def test_not_multiple_of_three():
    seq = "ATG" + "GCT" * 40 + "TAA" + "A"
    assert check_cds(seq) == (None, "not_multiple_of_three")


def test_ambiguous_bases():
    seq = "ATG" + "GCT" * 39 + "GNT" + "TAA"
    assert check_cds(seq) == (None, "ambiguous_bases")


def test_no_stop_codon():
    seq = "ATG" + "GCT" * 40
    assert check_cds(seq, require_stop=True) == (None, "no_stop_codon")
    assert check_cds(seq, require_stop=False)[1] == "ok"


def test_internal_stop():
    seq = "ATG" + "GCT" * 20 + "TGA" + "GCT" * 19 + "TAA"
    assert check_cds(seq) == (None, "internal_stop")


def test_bad_start():
    seq = "CCC" + "GCT" * 40 + "TAA"
    assert check_cds(seq, require_start=True) == (None, "bad_start")
    assert check_cds(seq, require_start=False)[1] == "ok"


def test_every_reason_is_declared():
    produced = set()
    cases = ["ATGGCTTAA", "ATG" + "GCT" * 40 + "TAAA",
             "ATG" + "GCT" * 39 + "GNT" + "TAA", "ATG" + "GCT" * 40,
             "ATG" + "GCT" * 20 + "TGA" + "GCT" * 19 + "TAA",
             "CCC" + "GCT" * 40 + "TAA", "ATG" + "GCT" * 40 + "TAA"]
    for seq in cases:
        produced.add(check_cds(seq, require_start=True)[1])
    assert produced == set(REASONS)


def test_rna_and_lowercase_are_accepted():
    seq = ("aug" + "gcu" * 40 + "uaa")
    cleaned, reason = check_cds(seq)
    assert reason == "ok"
    assert cleaned.startswith("ATGGCT")


def test_clean_cds_wrapper_returns_sequence_only():
    seq = "ATG" + "GCT" * 40 + "TAA"
    assert clean_cds(seq) == check_cds(seq)[0]
