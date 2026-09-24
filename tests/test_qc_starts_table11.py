"""Start-codon handling in check_cds, NCBI translation table 11.

Added 2026-09-24 with the fix that made STARTS_TABLE11 the default. Table 11
(Bacterial, Archaeal and Plant Plastid) marks seven codons as initiators:
ATG, GTG, TTG, ATT, ATC, ATA and CTG. The previous default accepted only the
first three, so a CDS beginning ATA was reported as 'bad_start' even though
AMRFinderPlus had already called it a 100 per cent protein-level allele match.
"""
import pytest

from codonamr.genetic_code import STARTS_CONSERVATIVE, STARTS_TABLE11
from codonamr.qc import check_cds

BODY = "AAA" * 40           # 40 lysine codons, no internal stop
STOP = "TAA"


def cds(start):
    return start + BODY + STOP


@pytest.mark.parametrize("start", STARTS_TABLE11)
def test_all_table11_initiators_accepted_by_default(start):
    cleaned, reason = check_cds(cds(start))
    assert reason == "ok"
    assert cleaned == start + BODY          # terminal stop removed


@pytest.mark.parametrize("start", ("ATT", "ATC", "ATA", "CTG"))
def test_the_four_new_initiators_were_rejected_before(start):
    assert check_cds(cds(start), starts=STARTS_CONSERVATIVE)[1] == "bad_start"
    assert check_cds(cds(start))[1] == "ok"


def test_a_genuine_non_initiator_is_still_rejected():
    # AAA is not an initiation codon in any translation table
    assert check_cds(cds("AAA"))[1] == "bad_start"


def test_conservative_set_still_reproducible():
    for start in ("ATG", "GTG", "TTG"):
        assert check_cds(cds(start), starts=STARTS_CONSERVATIVE)[1] == "ok"


def test_start_change_cannot_turn_ok_into_a_rejection():
    # anything the old default accepted, the new default must also accept
    for start in STARTS_CONSERVATIVE:
        assert check_cds(cds(start))[1] == "ok"


def test_require_start_false_ignores_the_set_entirely():
    assert check_cds(cds("AAA"), require_start=False)[1] == "ok"
