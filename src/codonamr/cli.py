"""Command line interface.

The headline command is

    codonamr run --ids ids.txt --out results/

and it is meant to need nothing else. The one thing it cannot supply for the
user is a contact address, because NCBI requires one with every request, so
that comes from ``--email`` or from the ``CODONAMR_EMAIL`` environment
variable. Set it once in a shell profile and the command above is the whole
interface.

Subcommands
-----------
run       fetch, QC, metrics, host comparison, controls, randomisation, report
metrics   per-sequence indices only, no network and no host comparison
validate  run the control validation and print the table
version   print the version
"""
import argparse
import os
import sys
from pathlib import Path

__all__ = ["main", "build_parser"]

VERSION = "0.1.0"

_EMAIL_HELP = (
    "NCBI requires a contact address with every E-utilities request. Pass "
    "--email you@example.org, or set CODONAMR_EMAIL once in your shell "
    "profile. Local FASTA input with --fasta needs no email."
)


def _read_ids(path_or_ids):
    if path_or_ids is None:
        return None
    out = []
    for item in path_or_ids:
        p = Path(item)
        if p.exists():
            out.extend(w for w in p.read_text().split() if not w.startswith("#"))
        else:
            out.extend(x for x in item.replace(",", " ").split() if x)
    return out


def build_parser():
    p = argparse.ArgumentParser(
        prog="codonamr",
        description="Codon usage analysis of resistance genes, with the "
                    "controls that decide whether the result means anything.")
    p.add_argument("--version", action="version", version="codonamr " + VERSION)
    sub = p.add_subparsers(dest="command", metavar="{run,metrics,validate,version}")

    r = sub.add_parser("run", help="full analysis, identifiers in and a report out")
    r.add_argument("--ids", nargs="+", metavar="FILE_OR_ID",
                   help="file with one accession per line, or accessions")
    r.add_argument("--fasta", nargs="+", metavar="FILE",
                   help="local CDS FASTA, instead of or alongside --ids")
    r.add_argument("--out", "-o", default="codonamr_results", metavar="DIR",
                   help="output directory (default: %(default)s)")
    r.add_argument("--email", default=None, help=_EMAIL_HELP)
    r.add_argument("--api-key", default=None,
                   help="NCBI API key, raising the rate limit from 3/s to 10/s")
    r.add_argument("--host", default=None, metavar="TAXID_OR_NAME_OR_FASTA",
                   help="override the host inferred from the records")
    r.add_argument("--db", default="auto", choices=("auto", "nuccore", "protein"),
                   help="which NCBI database the accessions belong to")
    r.add_argument("--no-folding", action="store_true",
                   help="skip mRNA folding even if ViennaRNA is installed")
    r.add_argument("--permutations", type=int, default=100, metavar="N",
                   help="permutations per gene (default: %(default)s)")
    r.add_argument("--folding-permutations", type=int, default=20, metavar="N",
                   help="how many of those permutations are also folded; "
                        "folding is the expensive step (default: %(default)s)")
    r.add_argument("--length-matched", type=int, default=100, metavar="N",
                   help="host genes drawn per length-matched null")
    r.add_argument("--seed", type=int, default=0,
                   help="master seed; the same seed gives the same output")
    r.add_argument("--threads", type=int, default=None, metavar="N",
                   help="worker processes; 1 runs everything in this process")
    r.add_argument("--cache-dir", default=None,
                   help="where downloads are cached (default: $CODONAMR_CACHE, "
                        "else the per-user cache directory)")
    r.add_argument("--trnascan", default=None, metavar="FILE",
                   help="tRNAscan-SE output for the host, enabling tAI")
    r.add_argument("--temperature", type=float, default=None, metavar="C",
                   help="folding temperature in Celsius (ViennaRNA default 37)")
    r.add_argument("--fold-window", type=int, default=40)
    r.add_argument("--fold-step", type=int, default=3)
    r.add_argument("--min-len", type=int, default=90, metavar="NT",
                   help="shortest acceptable CDS (default: %(default)s nt)")
    r.add_argument("--require-start", action="store_true",
                   help="reject sequences that do not begin at a start codon")
    r.add_argument("--no-require-stop", action="store_true",
                   help="accept sequences with no terminal stop codon")
    r.add_argument("--host-max-cds", type=int, default=None, metavar="N",
                   help="use only the first N host CDS, for a quick look")
    r.add_argument("--no-validate", action="store_true",
                   help="skip the control validation step")
    r.add_argument("--quiet", "-q", action="store_true")

    m = sub.add_parser("metrics", help="per-sequence indices only, offline")
    m.add_argument("--fasta", nargs="+", required=False, metavar="FILE")
    m.add_argument("--ids", nargs="+", metavar="FILE_OR_ID")
    m.add_argument("--email", default=None, help=_EMAIL_HELP)
    m.add_argument("--api-key", default=None)
    m.add_argument("--cache-dir", default=None)
    m.add_argument("--out", "-o", default="-", metavar="FILE",
                   help="TSV to write, or - for standard output")
    m.add_argument("--min-len", type=int, default=90)
    m.add_argument("--require-start", action="store_true")
    m.add_argument("--no-require-stop", action="store_true")
    m.add_argument("--reference", default=None, metavar="FASTA",
                   help="reference-set FASTA for CAI; without it CAI is NA")

    v = sub.add_parser("validate",
                       help="check that the controls recover what they should")
    v.add_argument("--n", type=int, default=200, help="genes per replicate")
    v.add_argument("--effect", type=float, default=0.3,
                   help="true partial correlation to inject")
    v.add_argument("--confound", type=float, default=0.8,
                   help="strength of the composition confounder")
    v.add_argument("--skew", type=float, default=0.5,
                   help="fraction of the confounder carried by AT skew, which "
                        "a GC-only control cannot see")
    v.add_argument("--replicates", type=int, default=500)
    v.add_argument("--alpha", type=float, default=0.05)
    v.add_argument("--seed", type=int, default=0)
    v.add_argument("--out", "-o", default=None, metavar="FILE",
                   help="also write the table as TSV")

    sub.add_parser("version", help="print the version")
    return p


def _cmd_run(args):
    from .pipeline import analyse
    email = args.email or os.environ.get("CODONAMR_EMAIL")
    if args.ids and not email:
        print("codonamr: --ids needs an email address. " + _EMAIL_HELP,
              file=sys.stderr)
        return 2
    if not args.ids and not args.fasta:
        print("codonamr: nothing to do, pass --ids or --fasta", file=sys.stderr)
        return 2
    result = analyse(
        ids=_read_ids(args.ids), fasta=args.fasta, out_dir=args.out,
        email=email, api_key=args.api_key or os.environ.get("CODONAMR_API_KEY"),
        host=args.host, folding=not args.no_folding,
        n_permutations=args.permutations, seed=args.seed, threads=args.threads,
        cache_dir=args.cache_dir, db=args.db, min_len=args.min_len,
        require_stop=not args.no_require_stop, require_start=args.require_start,
        trnascan=args.trnascan, temperature=args.temperature,
        fold_window=args.fold_window, fold_step=args.fold_step,
        n_folding=args.folding_permutations,
        n_length_matched=args.length_matched, host_max_cds=args.host_max_cds,
        validate=not args.no_validate, verbose=not args.quiet)
    m = result["manifest"]
    print("codonamr: %d/%d sequences analysed, host %s, report at %s"
          % (m["n_passed_qc"], m["n_requested"], m["host"] or "none",
             Path(args.out) / "report.md"))
    for n in result["notes"]:
        print("  note: %s" % n)
    return 0


def _cmd_metrics(args):
    import pandas as pd

    from .fetch import fetch_cds, read_fasta
    from .metrics import build_w, cai, enc, fop, gc3s, gc_content, scuo
    from .qc import check_cds

    records = []
    for p in (args.fasta or []):
        records.extend(read_fasta(p))
    if args.ids:
        email = args.email or os.environ.get("CODONAMR_EMAIL")
        if not email:
            print("codonamr: --ids needs an email address. " + _EMAIL_HELP,
                  file=sys.stderr)
            return 2
        records.extend(fetch_cds(_read_ids(args.ids), email=email,
                                 api_key=args.api_key,
                                 cache_dir=args.cache_dir, db="auto"))
    if not records:
        print("codonamr: nothing to do, pass --fasta or --ids", file=sys.stderr)
        return 2

    w = opt = None
    if args.reference:
        ref = [s for s in (check_cds(r.sequence, min_len=args.min_len,
                                     require_stop=False, require_start=False)[0]
                           for r in read_fasta(args.reference)) if s]
        if ref:
            from .host import optimal_codons
            from .metrics import rscu
            w = build_w(ref)
            opt = optimal_codons(rscu("".join(ref)))

    rows = []
    for r in records:
        seq, reason = check_cds(r.sequence, min_len=args.min_len,
                                require_stop=not args.no_require_stop,
                                require_start=args.require_start)
        row = {"id": r.id, "qc": reason, "organism": r.organism,
               "product": r.product}
        if seq:
            row.update({"n_codons": len(seq) // 3, "gc": gc_content(seq),
                        "gc3s": gc3s(seq), "enc": enc(seq), "scuo": scuo(seq),
                        "cai": cai(seq, w) if w else None,
                        "fop": fop(seq, opt) if opt else None})
        rows.append(row)
    df = pd.DataFrame(rows)
    if args.out in ("-", None):
        print(df.to_csv(sep="\t", index=False, float_format="%.6g",
                        na_rep="NA"), end="")
    else:
        from .pipeline import write_tsv
        write_tsv(df, args.out)
        print("codonamr: wrote %s" % args.out)
    return 0


def _cmd_validate(args):
    from .controls import false_positive_rate, injection_recovery
    import pandas as pd

    rec = injection_recovery(n=args.n, effect=args.effect,
                             confound=args.confound, skew=args.skew,
                             n_replicates=args.replicates, seed=args.seed,
                             alpha=args.alpha)
    rec.insert(0, "experiment", "injection_recovery")
    fpr = false_positive_rate(n=args.n, confound=args.confound, skew=args.skew,
                              n_replicates=args.replicates, seed=args.seed + 1,
                              alpha=args.alpha)
    fpr.insert(0, "experiment", "false_positive_rate")
    table = pd.concat([rec, fpr], ignore_index=True)
    print()
    print(table.to_string(index=False, float_format=lambda v: "%.4f" % v))
    print()
    print("injection_recovery: mean_r should be close to true_partial_r for a "
          "control that leaves real signal intact.")
    print("false_positive_rate: rejection_rate should be close to alpha = %.2f "
          "for a control that removes the confounding." % args.alpha)
    if args.out:
        from .pipeline import write_tsv
        write_tsv(table, args.out)
        print("wrote %s" % args.out)
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command in (None, "version"):
        if args.command is None:
            parser.print_help()
            return 0
        print("codonamr " + VERSION)
        return 0
    if args.command == "run":
        return _cmd_run(args)
    if args.command == "metrics":
        return _cmd_metrics(args)
    if args.command == "validate":
        return _cmd_validate(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
