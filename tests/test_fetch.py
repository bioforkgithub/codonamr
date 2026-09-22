"""Sequence input. Everything here is offline except the tests marked network."""
import os

import pytest

from codonamr.fetch import (CdsRecord, default_cache_dir, guess_db,
                            parse_cds_fasta, read_fasta)

pytestmark = []


def test_guess_db_recognises_refseq_prefixes():
    assert guess_db("WP_000027057.1") == "protein"
    assert guess_db("NP_414542.1") == "protein"
    assert guess_db("NC_000913.3") == "nuccore"
    assert guess_db("NG_048025.1") == "nuccore"
    assert guess_db("AAA12345.1") == "protein"
    assert guess_db("J01749") == "nuccore"


def test_read_fasta_parses_bracketed_tags(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text(">gene1 [organism=Escherichia coli] [gene=blaTEM] "
                 "[protein=beta-lactamase]\nATGGCTTAA\n"
                 ">gene2 plain description\nATGGGGTAA\n")
    recs = read_fasta(p)
    assert [r.id for r in recs] == ["gene1", "gene2"]
    assert recs[0].organism == "Escherichia coli"
    assert recs[0].gene == "blaTEM"
    assert recs[0].product == "beta-lactamase"
    assert recs[1].organism is None
    assert recs[0].sequence == "ATGGCTTAA"
    assert recs[0].source == "local"


def test_read_fasta_organism_override(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text(">a\nATGTAA\n")
    recs = read_fasta(p, organism="Klebsiella pneumoniae", taxid="573")
    assert recs[0].taxid == "573"


def test_parse_cds_fasta_reads_ncbi_cds_headers(tmp_path):
    p = tmp_path / "cds.fna"
    p.write_text(
        ">lcl|NC_000913.3_cds_NP_414542.1_1 [gene=thrA] [locus_tag=b0002] "
        "[protein=bifunctional aspartokinase] [protein_id=NP_414542.1] "
        "[location=337..2799] [gbkey=CDS]\nATGGCTTAA\n"
        ">lcl|NC_000913.3_cds_2 [gene=ghost] [pseudo=true] [gbkey=CDS]\n"
        "ATGGGGTAA\n")
    with open(p) as fh:
        recs = parse_cds_fasta(fh)
    assert len(recs) == 1
    assert recs[0].gene == "thrA"
    assert recs[0].product == "bifunctional aspartokinase"
    assert recs[0].extra["locus_tag"] == "b0002"


def test_cds_record_round_trips_to_a_dict():
    r = CdsRecord(id="x", sequence="ATG")
    assert r.as_dict()["id"] == "x"
    assert r.as_dict()["sequence"] == "ATG"


def test_default_cache_dir_honours_the_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("CODONAMR_CACHE", str(tmp_path / "cache"))
    assert default_cache_dir() == tmp_path / "cache"
    monkeypatch.delenv("CODONAMR_CACHE")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    assert default_cache_dir() == tmp_path / "codonamr"


def test_fetch_cds_requires_an_email():
    from codonamr.fetch import fetch_cds
    with pytest.raises(ValueError, match="email"):
        fetch_cds(["NG_048025.1"], email=None)


@pytest.mark.network
def test_fetch_cds_from_ncbi(tmp_path):
    """Needs the network. Run with `pytest -m network`."""
    from codonamr.fetch import fetch_cds
    email = os.environ.get("CODONAMR_EMAIL")
    if not email:
        pytest.skip("set CODONAMR_EMAIL to run the network tests")
    recs = fetch_cds(["NG_048025.1"], email=email, cache_dir=tmp_path)
    assert recs
    assert recs[0].sequence
    assert recs[0].organism
    # second call must be served from the cache, so it must still work offline
    again = fetch_cds(["NG_048025.1"], email=email, cache_dir=tmp_path)
    assert [r.sequence for r in again] == [r.sequence for r in recs]


@pytest.mark.network
def test_resolve_host_genome(tmp_path):
    from codonamr.fetch import resolve_host_genome
    email = os.environ.get("CODONAMR_EMAIL")
    if not email:
        pytest.skip("set CODONAMR_EMAIL to run the network tests")
    got = resolve_host_genome("562", email=email, cache_dir=tmp_path,
                              max_cds=50)
    assert got and got["cds"]
    assert "Escherichia" in got["assembly"].organism
