"""Host codon tables, reference-set selection, tAI and gene-to-host distance."""
import pytest

from codonamr.genetic_code import CODON2AA, FAMILY
from codonamr.host import (ILE2_KEY, cufs, host_codon_table, optimal_codons,
                           reference_set, rscu_distance, tai, tai_weights,
                           trna_copy_numbers_from_tRNAscan)
from codonamr.metrics import rscu

from conftest import one_codon_per_aa, random_cds, uniform_usage


def test_host_codon_table_pools_counts(host_cds):
    table = host_codon_table(host_cds)
    assert table["n_genes"] == len(host_cds)
    assert table["n_codons"] == sum(table["counts"].values())
    assert set(table["w"]) >= {c for c in CODON2AA if CODON2AA[c] != "*"}


def test_optimal_codons_are_those_above_rscu_one():
    table = rscu(one_codon_per_aa(repeats=10))
    opt = optimal_codons(table)
    for aa, fam in FAMILY.items():
        if len(fam) == 1:
            assert fam[0] not in opt
        else:
            assert fam[0] in opt
            assert all(c not in opt for c in fam[1:])


def test_reference_set_keeps_ribosomal_proteins(rng):
    seqs = [random_cds(rng, 200) for _ in range(30)]
    products = ["50S ribosomal protein L%d" % i for i in range(30)]
    assert len(reference_set(seqs, products, min_genes=5)) == 30


def test_reference_set_excludes_ribosomal_protein_modifying_enzymes(rng):
    """prmA, rimK, rimO, roxA and ycaO name a ribosomal protein but are enzymes."""
    decoys = [
        ("ribosomal protein L11 methyltransferase", "prmA"),
        ("50S ribosomal protein L3 N(5)-glutamine methyltransferase", "prmB"),
        ("ribosomal protein S6--L-glutamate ligase", "rimK"),
        ("ribosomal protein S12 methylthiotransferase", "rimO"),
        ("50S ribosomal protein L16 3-hydroxylase", "roxA"),
        ("ribosomal protein S12 methylthiotransferase accessory factor", "ycaO"),
        ("ribosomal protein alanine N-acetyltransferase", "rimI"),
        ("16S rRNA (guanine(527)-N(7))-methyltransferase", "rsmG"),
        ("ribosome maturation factor RimP", "rimP"),
        ("ribosome-binding ATPase YchF", "ychF"),
    ]
    real = [("30S ribosomal protein S%d" % i, "rpsA") for i in range(25)]
    seqs = [random_cds(rng, 200) for _ in range(len(decoys) + len(real))]
    products = [p for p, _ in decoys + real]
    genes = [g for _, g in decoys + real]
    kept = reference_set(seqs, products, genes, min_genes=5)
    assert len(kept) == len(real)
    assert all(s in seqs[len(decoys):] for s in kept)


def test_reference_set_returns_nothing_when_too_few_hits(rng):
    seqs = [random_cds(rng, 200) for _ in range(5)]
    products = ["50S ribosomal protein L%d" % i for i in range(5)]
    assert reference_set(seqs, products, min_genes=20) == []


def test_tai_weights_are_normalised_and_bounded():
    counts = {"AGC": 4, "GCT": 2, "CAT": 3, "GAT": 2, "TAC": 1, "TTC": 5,
              ILE2_KEY: 1}
    w = tai_weights(counts)
    assert len(w) == 61
    assert max(w.values()) == pytest.approx(1.0)
    assert min(w.values()) > 0


def test_tai_of_a_gene_using_only_the_best_decoded_codon_is_one():
    counts = {"AGC": 10, ILE2_KEY: 1}
    w = tai_weights(counts)
    best = max((v, c) for c, v in w.items())[1]
    seq = best * 50
    assert tai(seq, counts) == pytest.approx(1.0, abs=1e-12)


def test_tai_degrades_gracefully_without_trna_data():
    assert tai("ATGGCTGCTGCT", {}) is None
    assert tai("ATGGCTGCTGCT", None) is None


def test_tai_wobble_penalties_match_the_published_constants():
    """Each wobble pairing costs exactly the s value of dos Reis et al. (2004)."""
    # anticodon AGC is A34, Watson-Crick on GCT and inosine wobble on GCC
    w = tai_weights({"AGC": 10})
    assert w["GCT"] > w["GCC"]
    assert w["GCC"] / w["GCT"] == pytest.approx(1 - 0.28, abs=1e-9)
    # anticodon GGC is G34, Watson-Crick on GCC and G:U wobble on GCT
    w = tai_weights({"GGC": 10})
    assert w["GCC"] > w["GCT"]
    assert w["GCT"] / w["GCC"] == pytest.approx(1 - 0.41, abs=1e-9)
    # anticodon TGC is U34, Watson-Crick on GCA and U:G wobble on GCG
    w = tai_weights({"TGC": 10})
    assert w["GCG"] / w["GCA"] == pytest.approx(1 - 0.68, abs=1e-9)


def test_trnascan_parser(tmp_path):
    text = """Sequence\t\ttRNA \tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t
Name    \ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote
--------\t------\t-----\t----\t----\t-----\t-----\t----\t------\t------
chr1\t1\t100\t176\tAla\tTGC\t0\t0\t85.1\t
chr1\t2\t300\t376\tAla\tTGC\t0\t0\t80.0\t
chr1\t3\t500\t576\tIle\tCAT\t0\t0\t70.0\t
chr1\t4\t700\t776\tMet\tCAT\t0\t0\t90.0\t
chr1\t5\t900\t976\tUndet\tNNN\t0\t0\t20.0\t
chr1\t6\t990\t1066\tLeu\tCAA\t0\t0\t55.0\tpseudo
"""
    p = tmp_path / "trna.out"
    p.write_text(text)
    counts = trna_copy_numbers_from_tRNAscan(p)
    assert counts["TGC"] == 2
    assert counts["CAT"] == 1          # Met only
    assert counts[ILE2_KEY] == 1       # lysidine tRNA-Ile kept separate
    assert "NNN" not in counts
    assert "CAA" not in counts


def test_cufs_is_zero_for_identical_usage_and_positive_otherwise():
    a = rscu(one_codon_per_aa(repeats=10))
    assert cufs(a, a) == pytest.approx(0.0, abs=1e-12)
    b = rscu(uniform_usage(repeats=10))
    assert 0.0 < cufs(a, b) <= 1.0


def test_cufs_is_symmetric(rng):
    a = rscu(random_cds(rng, 400, bias=0.8))
    b = rscu(random_cds(rng, 400, bias=0.2))
    assert cufs(a, b) == pytest.approx(cufs(b, a), abs=1e-12)


def test_rscu_distance_handles_missing_families():
    a = rscu("ATGGCTGCCTAA")      # only Ala is informative
    b = rscu(uniform_usage(repeats=5))
    assert rscu_distance(a, a) == pytest.approx(0.0)
    assert rscu_distance(a, b) is not None
    assert rscu_distance(rscu("ATGATG"), b) is None


def test_rscu_distance_metrics_agree_on_identity(rng):
    a = rscu(random_cds(rng, 600, bias=0.7, stop=""))
    for metric in ("euclidean", "manhattan", "cosine", "correlation"):
        assert rscu_distance(a, a, metric) == pytest.approx(0.0, abs=1e-9)
    with pytest.raises(ValueError):
        rscu_distance(a, a, "nonsense")


def test_rscu_correlation_distance_is_undefined_for_a_flat_table():
    """Every RSCU equal to 1 has no variance, so a correlation is meaningless."""
    flat = rscu(uniform_usage(repeats=5))
    assert rscu_distance(flat, flat, "correlation") is None
    assert rscu_distance(flat, flat, "euclidean") == pytest.approx(0.0)
