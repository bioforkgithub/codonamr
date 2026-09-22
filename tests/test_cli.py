"""The command line, including the headline one-liner."""
import pytest

from codonamr.cli import build_parser, main


def test_help_exits_cleanly(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert "run" in capsys.readouterr().out


def test_version(capsys):
    assert main(["version"]) == 0
    assert "0.1.0" in capsys.readouterr().out


def test_no_arguments_prints_help(capsys):
    assert main([]) == 0
    assert "usage" in capsys.readouterr().out


def test_run_without_input_is_an_error(capsys):
    assert main(["run"]) == 2
    assert "pass --ids or --fasta" in capsys.readouterr().err


def test_ids_without_an_email_explains_itself(capsys, monkeypatch):
    monkeypatch.delenv("CODONAMR_EMAIL", raising=False)
    assert main(["run", "--ids", "NG_048025.1"]) == 2
    assert "email" in capsys.readouterr().err.lower()


def test_metrics_subcommand_on_a_local_fasta(tmp_path, capsys):
    p = tmp_path / "in.fasta"
    p.write_text(">g1 [protein=test]\n" + "ATG" + "GCT" * 40 + "TAA" + "\n")
    assert main(["metrics", "--fasta", str(p)]) == 0
    out = capsys.readouterr().out
    assert "enc" in out and "g1" in out
    assert "\t" in out


def test_metrics_subcommand_writes_a_file(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text(">g1\n" + "ATG" + "GCT" * 40 + "TAA" + "\n")
    out = tmp_path / "m.tsv"
    assert main(["metrics", "--fasta", str(p), "--out", str(out)]) == 0
    assert "enc" in out.read_text()


def test_validate_subcommand_prints_the_table(capsys):
    assert main(["validate", "--n", "80", "--replicates", "20"]) == 0
    text = capsys.readouterr().out
    assert "injection_recovery" in text
    assert "false_positive_rate" in text
    assert "per-nucleotide composition" in text


def test_run_subcommand_end_to_end(tmp_path):
    from conftest import random_cds
    import random
    rng = random.Random(1)
    host = tmp_path / "host.fasta"
    with open(host, "w") as fh:
        for i in range(40):
            fh.write(">h%d [protein=30S ribosomal protein S%d] [gene=rpsA]\n%s\n"
                     % (i, i, random_cds(rng, 200, bias=0.9)))
        for i in range(160):
            fh.write(">o%d [protein=hypothetical protein]\n%s\n"
                     % (i, random_cds(rng, rng.randrange(80, 400), bias=0.5)))
    query = tmp_path / "q.fasta"
    with open(query, "w") as fh:
        for i in range(3):
            fh.write(">q%d [protein=beta-lactamase]\n%s\n"
                     % (i, random_cds(rng, 250, bias=0.7)))
    out = tmp_path / "res"
    code = main(["run", "--fasta", str(query), "--host", str(host),
                 "--out", str(out), "--no-folding", "--permutations", "5",
                 "--length-matched", "10", "--no-validate", "--threads", "1",
                 "--quiet"])
    assert code == 0
    assert (out / "report.md").exists()
    assert (out / "metrics.tsv").exists()


def test_parser_declares_all_four_subcommands():
    parser = build_parser()
    actions = [a for a in parser._actions if hasattr(a, "choices") and a.choices]
    names = set()
    for a in actions:
        names |= set(a.choices)
    assert {"run", "metrics", "validate", "version"} <= names
