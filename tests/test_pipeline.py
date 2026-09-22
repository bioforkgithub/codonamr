"""End to end, offline, on synthetic data.

The pipeline is tested with local FASTA input and a local host, so the whole
suite runs with no network. What is asserted is that every promised file
appears, that the same seed gives the same numbers, that rejected sequences are
reported rather than dropped, and that a single-process run and a multi-process
run agree.
"""
import json
import random
from pathlib import Path

import pandas as pd
import pytest

from codonamr.pipeline import analyse

from conftest import random_cds


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, desc, seq in records:
            fh.write(">%s %s\n%s\n" % (name, desc, seq))
    return path


@pytest.fixture
def host_fasta(tmp_path):
    rng = random.Random(4242)
    recs = []
    for i in range(30):
        recs.append(("host_rp_%d" % i,
                     "[organism=Testus bacterius] [gene=rpsA] "
                     "[protein=30S ribosomal protein S%d]" % i,
                     random_cds(rng, rng.randrange(90, 260), bias=0.9)))
    for i in range(6):
        recs.append(("host_decoy_%d" % i,
                     "[organism=Testus bacterius] [gene=prmA] "
                     "[protein=ribosomal protein L11 methyltransferase]",
                     random_cds(rng, 250, bias=0.2)))
    for i in range(240):
        recs.append(("host_gene_%d" % i,
                     "[organism=Testus bacterius] [gene=g%d] "
                     "[protein=hypothetical protein]" % i,
                     random_cds(rng, rng.randrange(60, 600), bias=0.55)))
    return _write_fasta(tmp_path / "host.fasta", recs)


@pytest.fixture
def query_fasta(tmp_path):
    """Twelve usable genes and one of each QC failure mode.

    Twelve because the composition control needs at least eight genes with
    complete values before it will report anything, and a fixture that never
    exercises it would leave the most important table untested.
    """
    rng = random.Random(7)
    recs = []
    for i in range(12):
        recs.append(("amr_%02d" % (i + 1),
                     "[organism=Testus bacterius] [protein=beta-lactamase %d]" % i,
                     random_cds(rng, rng.randrange(140, 420),
                                bias=0.3 + 0.05 * i)))
    recs += [
        ("bad_frame", "[organism=Testus bacterius] [protein=broken]",
         random_cds(rng, 200) + "A"),
        ("bad_short", "[organism=Testus bacterius] [protein=too short]",
         "ATGGCTTAA"),
        ("bad_ambiguous", "[organism=Testus bacterius] [protein=ambiguous]",
         "ATG" + "GCT" * 50 + "GNT" + "TAA"),
    ]
    return _write_fasta(tmp_path / "query.fasta", recs)


@pytest.fixture
def run(tmp_path, query_fasta, host_fasta):
    return analyse(fasta=query_fasta, host=host_fasta,
                   out_dir=tmp_path / "out", folding=False,
                   n_permutations=15, n_length_matched=30, seed=0, threads=1,
                   verbose=False)


EXPECTED_FILES = [
    "qc.tsv", "qc_summary.tsv", "metrics.tsv", "rscu_per_gene.tsv",
    "aa_composition.tsv", "codon_counts.tsv", "host_codon_table.tsv",
    "controls_length_matched.tsv", "controls_composition.tsv",
    "randomisation.tsv", "per_gene_all.tsv", "validation.tsv",
    "manifest.json", "report.md",
]


def test_every_intermediate_is_written_as_tsv(run):
    out = Path(run["out_dir"])
    for name in EXPECTED_FILES:
        assert (out / name).exists(), name
    figures = list((out / "figures").glob("*.png"))
    assert len(figures) >= 5
    assert len(list((out / "figures").glob("*.pdf"))) >= 5


def test_qc_reports_every_rejection_with_a_reason(run):
    qc = run["qc"]
    reasons = dict(zip(qc["id"], qc["reason"]))
    assert reasons["bad_frame"] == "not_multiple_of_three"
    assert reasons["bad_short"] == "too_short"
    assert reasons["bad_ambiguous"] == "ambiguous_bases"
    assert reasons["amr_01"] == "ok"
    assert reasons["amr_12"] == "ok"
    assert run["manifest"]["n_requested"] == 15
    assert run["manifest"]["n_passed_qc"] == 12
    assert int(run["qc_summary"].set_index("reason")
               .loc["ok", "n_sequences"]) == 12


def test_report_lists_a_metric_that_could_not_be_computed(run):
    """tAI has no tRNA data here, so it must appear in the missing table."""
    from codonamr.report import _missing_table
    miss = _missing_table(run["per_gene"]).set_index("metric")
    assert "tai" in miss.index
    assert int(miss.loc["tai", "n_missing"]) == 12
    assert "tRNA" in miss.loc["tai", "reason"]
    text = (Path(run["out_dir"]) / "report.md").read_text()
    assert "no tRNA gene copy numbers supplied" in text


def test_report_states_the_rejection_counts(run):
    text = (Path(run["out_dir"]) / "report.md").read_text()
    assert "Quality control" in text
    assert "not_multiple_of_three" in text
    assert "Missing values" in text
    assert "Nothing was dropped silently" in text
    # the report must contain no em-dash, written here as an escape so
    # that this file does not contain one either
    assert "\u2014" not in text


def test_metrics_are_present_for_every_kept_sequence(run):
    m = run["metrics"]
    assert list(m["id"]) == ["amr_%02d" % i for i in range(1, 13)]
    for col in ("gc", "gc3s", "enc", "scuo", "cai_reference", "fop",
                "cufs_host", "n_codons"):
        assert m[col].notna().all(), col
    assert m["tai"].isna().all()      # no tRNA data was supplied


def test_host_reference_set_excluded_the_decoys(run):
    assert run["manifest"]["reference_set_size"] == 30
    assert run["manifest"]["host_n_cds"] > 200


def test_composition_control_is_reported(run):
    cc = run["controls_composition"]
    assert len(cc)
    assert {"x", "y", "r_raw", "p_raw", "r_partial", "p_partial"} <= set(cc.columns)
    assert (cc["k"] == 3).all()          # A, C and G, with T dropped
    assert cc["r_partial"].notna().any()


def test_length_matched_null_is_length_matched(run):
    cl = run["controls_length_matched"]
    assert set(cl["metric"]) == {"enc", "scuo", "gc3s", "cai_reference"}
    assert (cl["n_null"] > 0).all()
    assert cl["pool_size"].min() > 0


def test_same_seed_gives_the_same_output(tmp_path, query_fasta, host_fasta):
    kw = dict(fasta=query_fasta, host=host_fasta, folding=False,
              n_permutations=15, n_length_matched=30, threads=1, verbose=False)
    analyse(out_dir=tmp_path / "a", seed=11, **kw)
    analyse(out_dir=tmp_path / "b", seed=11, **kw)
    for name in ("metrics.tsv", "randomisation.tsv",
                 "controls_length_matched.tsv", "validation.tsv"):
        assert (tmp_path / "a" / name).read_text() == \
               (tmp_path / "b" / name).read_text(), name


def test_a_different_seed_changes_the_randomisation(tmp_path, query_fasta,
                                                    host_fasta):
    kw = dict(fasta=query_fasta, host=host_fasta, folding=False,
              n_permutations=15, n_length_matched=30, threads=1, verbose=False,
              validate=False)
    a = analyse(out_dir=tmp_path / "a", seed=1, **kw)
    b = analyse(out_dir=tmp_path / "b", seed=2, **kw)
    assert not a["controls_length_matched"].equals(b["controls_length_matched"])
    # the metrics themselves are deterministic and must not move
    pd.testing.assert_frame_equal(a["metrics"], b["metrics"])


def test_threads_do_not_change_the_result(tmp_path, query_fasta, host_fasta):
    kw = dict(fasta=query_fasta, host=host_fasta, folding=False,
              n_permutations=15, n_length_matched=30, seed=3, verbose=False,
              validate=False)
    one = analyse(out_dir=tmp_path / "one", threads=1, **kw)
    many = analyse(out_dir=tmp_path / "many", threads=2, **kw)
    pd.testing.assert_frame_equal(one["randomisation"], many["randomisation"])
    pd.testing.assert_frame_equal(one["controls_length_matched"],
                                  many["controls_length_matched"])


def test_randomisation_reports_the_dinucleotide_residual(run):
    rand = run["randomisation"]
    assert (rand["dinuc_residual_mean"] >= 0).all()
    assert rand["dinuc_residual_zero_fraction"].between(0, 1).all()
    assert rand["junction_chi2_p"].between(0, 1).all()
    assert rand["n_permutations"].eq(15).all()


def test_manifest_records_what_was_missing(run):
    manifest = json.loads((Path(run["out_dir"]) / "manifest.json").read_text())
    assert manifest["codonamr_version"] == "0.1.0"
    assert manifest["tai_available"] is False
    assert any("tAI" in n for n in manifest["notes"])


def test_analyse_refuses_an_empty_input(tmp_path):
    with pytest.raises(ValueError, match="no sequences"):
        analyse(out_dir=tmp_path / "x", verbose=False)


def test_analyse_reports_when_everything_fails_qc(tmp_path, host_fasta):
    bad = _write_fasta(tmp_path / "bad.fasta",
                       [("a", "", "ATGGCTTAA"), ("b", "", "ATGGGGTAA")])
    with pytest.raises(ValueError, match="failed QC"):
        analyse(fasta=bad, host=host_fasta, out_dir=tmp_path / "y",
                folding=False, verbose=False)
    assert (tmp_path / "y" / "qc_summary.tsv").exists()


@pytest.mark.skipif(not __import__("codonamr.folding", fromlist=["available"])
                    .available(), reason="ViennaRNA not installed")
def test_folding_path_runs_when_viennarna_is_present(tmp_path, query_fasta,
                                                     host_fasta):
    res = analyse(fasta=query_fasta, host=host_fasta,
                  out_dir=tmp_path / "fold", folding=True, n_permutations=3,
                  n_folding=2, n_length_matched=10, fold_step=20, seed=0,
                  threads=1, validate=False, verbose=False)
    rand = res["randomisation"]
    assert rand["dg_mean"].notna().all()
    assert (rand["dg_mean"] < 0).all()      # folding free energies are negative
    assert rand["dg_p"].between(0, 1).all()
    assert res["manifest"]["folding"] is True
