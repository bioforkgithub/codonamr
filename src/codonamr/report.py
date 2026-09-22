"""Markdown report and figures.

The report is Markdown on purpose. It is readable in a terminal, diffable in
version control, and it is the same text whether it is read on disk or rendered,
which a PDF is not. Every figure is written twice, as a 300 dpi PNG for
inclusion and as a vector PDF for a manuscript, and every figure has the TSV it
was drawn from sitting next to it.

The report states what was rejected and what came back missing. A metric that
returned ``None`` is listed with its count and its reason, because a pipeline
that quietly drops the sequences it cannot handle reports the properties of the
sequences it happened to like.

References for the plotted relationships
----------------------------------------
ENC against GC3s    Wright F (1990) Gene 87:23-29.
Neutrality plot     Sueoka N (1988) Proc Natl Acad Sci USA 85:2653-2657.
PR2 bias plot       Sueoka N (1995) J Mol Evol 40:318-325.
"""
import datetime as _dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from .genetic_code import CODON2AA, CODONS

__all__ = ["write_report", "wright_expected_enc", "plot_enc_gc3s",
           "plot_neutrality", "plot_pr2", "plot_rscu_heatmap",
           "plot_host_compatibility", "PALETTE"]

#: Okabe and Ito's qualitative palette, which is distinguishable under
#: deuteranopia, protanopia and tritanopia and survives greyscale printing.
#: Okabe M, Ito K (2008) Color Universal Design.
PALETTE = {
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "green": "#009E73",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "ink": "#1A1A1A",
    "muted": "#6E6E6E",
    "grid": "#D9D9D9",
}
SERIES = [PALETTE["blue"], PALETTE["vermillion"], PALETTE["green"],
          PALETTE["orange"], PALETTE["purple"], PALETTE["sky"]]

#: diverging ramp for RSCU, two hues with a neutral grey midpoint at RSCU = 1,
#: which is the value that means "no bias" and so has to be the neutral point
RSCU_CMAP = LinearSegmentedColormap.from_list(
    "rscu", [PALETTE["blue"], "#C9D7E4", "#EDEDEA", "#F0CDB4",
             PALETTE["vermillion"]])

_RC = {
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.edgecolor": PALETTE["muted"],
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": PALETTE["grid"],
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "xtick.color": PALETTE["ink"],
    "ytick.color": PALETTE["ink"],
    "savefig.bbox": "tight",
}


def _style():
    return plt.rc_context(_RC)


def _save(fig, out_dir, name):
    """Write a figure as 300 dpi PNG and as vector PDF, and return the paths."""
    figs = out_dir / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    png, pdf = figs / (name + ".png"), figs / (name + ".pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def _finish(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", color=PALETTE["ink"])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# ------------------------------------------------------------------ curves
def wright_expected_enc(gc3s):
    """Wright's (1990) expected ENC under no selection, given GC3s.

    ``Nc = 2 + s + 29 / (s^2 + (1 - s)^2)`` with ``s`` the GC3s value. Points on
    the curve are consistent with mutational bias alone; points well below it
    need an additional explanation, usually translational selection. Points
    above the curve are not evidence of anything except estimation noise, which
    is common in short genes.
    """
    s = np.asarray(gc3s, dtype=float)
    return 2.0 + s + 29.0 / (s ** 2 + (1.0 - s) ** 2)


# ------------------------------------------------------------------- plots
def plot_enc_gc3s(genes, host=None, out_dir=None, name="enc_vs_gc3s"):
    """ENC against GC3s with the Wright expected curve.

    ``genes`` and ``host`` are DataFrames with ``gc3s`` and ``enc`` columns.
    """
    with _style():
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        s = np.linspace(0.01, 0.99, 400)
        ax.plot(s, wright_expected_enc(s), color=PALETTE["ink"], lw=1.6,
                zorder=2, label="Wright expectation, no selection")
        if host is not None and len(host):
            ax.scatter(host["gc3s"], host["enc"], s=9, alpha=0.30,
                       color=PALETTE["muted"], linewidths=0, zorder=1,
                       label="host genes (n=%d)" % len(host))
        g = genes.dropna(subset=["gc3s", "enc"])
        ax.scatter(g["gc3s"], g["enc"], s=34, color=PALETTE["blue"],
                   edgecolors="white", linewidths=0.8, zorder=3,
                   label="query genes (n=%d)" % len(g))
        ax.set_xlim(0, 1)
        ax.set_ylim(18, 63)
        _finish(ax, "Codon bias against synonymous GC content",
                "GC3s", "Effective number of codons")
        ax.legend(loc="lower center", fontsize=8, ncol=1)
        return _save(fig, out_dir, name) if out_dir else fig


def plot_neutrality(genes, out_dir=None, name="neutrality"):
    """GC12 against GC3, the Sueoka (1988) neutrality plot.

    The regression slope estimates the fraction of variation attributable to
    genome-wide mutational pressure. A slope near 1 means positions 1 and 2
    track position 3, so mutation dominates; a slope near 0 means they do not,
    which is the signature of selection acting on the coding positions.
    """
    g = genes.dropna(subset=["gc12", "gc3"])
    with _style():
        fig, ax = plt.subplots(figsize=(5.0, 4.0))
        ax.scatter(g["gc3"], g["gc12"], s=32, color=PALETTE["blue"],
                   edgecolors="white", linewidths=0.8, zorder=3,
                   label="query genes (n=%d)" % len(g))
        if len(g) >= 3 and g["gc3"].std() > 0:
            b, a = np.polyfit(g["gc3"], g["gc12"], 1)
            xs = np.linspace(g["gc3"].min(), g["gc3"].max(), 50)
            r = float(np.corrcoef(g["gc3"], g["gc12"])[0, 1])
            ax.plot(xs, a + b * xs, color=PALETTE["vermillion"], lw=2,
                    zorder=4, label="slope %.3f, r %.2f" % (b, r))
        ax.plot([0, 1], [0, 1], color=PALETTE["muted"], lw=1, ls="--",
                zorder=2, label="slope 1, pure mutational drift")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        _finish(ax, "Neutrality plot", "GC3", "GC12")
        ax.legend(loc="upper left", fontsize=8)
        return _save(fig, out_dir, name) if out_dir else fig


def plot_pr2(genes, out_dir=None, name="pr2_bias"):
    """Parity rule 2 plot at fourfold-degenerate third positions.

    Under no bias, mutation and selection affect both strands equally and both
    coordinates sit at 0.5, the centre of the plot. Distance from the centre
    measures strand asymmetry, and the direction says which of A/T and G/C is
    favoured. This is the figure that shows why a GC-only control is not enough:
    a gene can sit at GC3 = 0.5 and still be far off centre here.
    """
    g = genes.dropna(subset=["pr2_a3_at3", "pr2_g3_gc3"])
    with _style():
        fig, ax = plt.subplots(figsize=(4.6, 4.4))
        ax.axhline(0.5, color=PALETTE["muted"], lw=1, ls="--", zorder=2)
        ax.axvline(0.5, color=PALETTE["muted"], lw=1, ls="--", zorder=2)
        ax.scatter(g["pr2_a3_at3"], g["pr2_g3_gc3"], s=34,
                   color=PALETTE["blue"], edgecolors="white", linewidths=0.8,
                   zorder=3, label="query genes (n=%d)" % len(g))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        _finish(ax, "PR2 bias at fourfold sites",
                "A3 / (A3 + T3)", "G3 / (G3 + C3)")
        ax.annotate("no strand bias", xy=(0.5, 0.5), xytext=(0.56, 0.56),
                    fontsize=8, color=PALETTE["muted"])
        return _save(fig, out_dir, name) if out_dir else fig


def plot_rscu_heatmap(rscu_table, host_rscu=None, out_dir=None,
                      name="rscu_heatmap", max_genes=40):
    """RSCU per gene, codons grouped by amino acid, diverging about 1.

    RSCU = 1 means the codon is used exactly as often as even usage would
    predict, so it is the neutral point of the colour scale and is drawn grey.
    Blue is under-use and orange is over-use. Flat grey cells mean the family
    was never used in that gene, which is not the same as RSCU 0 and is common
    in short genes.

    The colour scale saturates at the 98th percentile of the values shown.
    RSCU has no upper bound and a six-codon family seen three times can reach 6,
    so scaling to the maximum would compress every real difference into the
    middle of the ramp.
    """
    order = [c for c in CODONS if CODON2AA[c] != "*"]
    order.sort(key=lambda c: (CODON2AA[c], c))
    df = rscu_table.set_index("id") if "id" in rscu_table else rscu_table
    df = df.reindex(columns=[c for c in order if c in df.columns])
    truncated = len(df) > max_genes
    if truncated:
        df = df.iloc[:max_genes]
    rows = list(df.index)
    data = df.to_numpy(dtype=float)
    if host_rscu is not None:
        hv = np.array([[host_rscu.get(c) if host_rscu.get(c) is not None
                        else np.nan for c in df.columns]], dtype=float)
        data = np.vstack([hv, data])
        rows = ["HOST"] + rows
    with _style():
        h = max(2.4, 0.22 * len(rows) + 1.6)
        fig, ax = plt.subplots(figsize=(min(14, 0.17 * data.shape[1] + 2.5), h))
        ax.grid(False)
        # The upper end is the 98th percentile rather than the maximum, so a
        # single extreme cell, which a short gene produces easily from small
        # counts, does not flatten the rest of the scale. Values above it
        # saturate, and the colourbar label says so.
        if np.isfinite(data).any():
            vmax = float(np.nanpercentile(data, 98))
        else:
            vmax = 2.0
        norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=max(vmax, 1.01))
        masked = np.ma.masked_invalid(data)
        cmap = RSCU_CMAP.copy()
        cmap.set_bad(PALETTE["grid"])
        im = ax.imshow(masked, aspect="auto", cmap=cmap, norm=norm)
        ax.set_xticks(range(data.shape[1]))
        ax.set_xticklabels(df.columns, rotation=90, fontsize=5.5)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows, fontsize=6.5)
        prev, bounds = None, []
        for i, c in enumerate(df.columns):
            if CODON2AA[c] != prev:
                bounds.append(i - 0.5)
                prev = CODON2AA[c]
        for b in bounds[1:]:
            ax.axvline(b, color="white", lw=1.2)
        cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("RSCU (1 = even use, top 2% saturated)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        title = "Relative synonymous codon usage, codons grouped by amino acid"
        if truncated:
            title += " (first %d genes)" % max_genes
        _finish(ax, title, "", "")
        return _save(fig, out_dir, name) if out_dir else fig


def plot_host_compatibility(genes, host=None, out_dir=None,
                            name="host_compatibility"):
    """How well each gene fits its host's codon preferences.

    Left: CAI against the host reference set, for the query genes and for the
    host's own genes, which is the distribution a gene has to be read against.
    A gene at the low end of the host distribution is using codons the host
    translates slowly. Right: CAI against CUFS distance from the host codon
    table, which separates a gene that is merely unbiased from one that is
    biased in a direction the host does not share.
    """
    g = genes
    with _style():
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))
        ax = axes[0]
        ax.grid(axis="x", visible=False)
        if host is not None and len(host.dropna(subset=["cai_reference"])):
            ax.hist(host["cai_reference"].dropna(), bins=40, density=True,
                    color=PALETTE["muted"], alpha=0.45,
                    label="host genes (n=%d)" % host["cai_reference"].notna().sum())
        vals = g["cai_reference"].dropna() if "cai_reference" in g else pd.Series(dtype=float)
        for i, v in enumerate(vals):
            ax.axvline(v, color=PALETTE["blue"], lw=1.4, alpha=0.85,
                       label="query genes (n=%d)" % len(vals) if i == 0 else None)
        _finish(ax, "CAI against the host reference set", "CAI", "density")
        ax.legend(loc="upper left", fontsize=8)

        ax = axes[1]
        if "cufs_host" in g and g["cufs_host"].notna().any():
            sub = g.dropna(subset=["cai_reference", "cufs_host"])
            ax.scatter(sub["cai_reference"], sub["cufs_host"], s=34,
                       color=PALETTE["blue"], edgecolors="white",
                       linewidths=0.8, zorder=3,
                       label="query genes (n=%d)" % len(sub))
            worst = sub.nlargest(min(3, len(sub)), "cufs_host")
            for _, row in worst.iterrows():
                ax.annotate(str(row["id"])[:18],
                            (row["cai_reference"], row["cufs_host"]),
                            fontsize=7, color=PALETTE["muted"],
                            xytext=(3, 3), textcoords="offset points")
            ax.legend(loc="upper right", fontsize=8)
        else:
            ax.text(0.5, 0.5, "no host codon table available",
                    ha="center", va="center", color=PALETTE["muted"])
        _finish(ax, "Fit against distance from host usage",
                "CAI", "CUFS distance from host (0 = identical)")
        fig.tight_layout()
        return _save(fig, out_dir, name) if out_dir else fig


# ------------------------------------------------------------------ report
def _md_table(df, floatfmt="%.4g", max_rows=None):
    if df is None or not len(df):
        return "_no rows_\n"
    d = df.head(max_rows) if max_rows else df
    cols = list(d.columns)
    def cell(v):
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "NA"
        if isinstance(v, float):
            return floatfmt % v
        return str(v)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, row in d.iterrows():
        lines.append("| " + " | ".join(cell(row[c]) for c in cols) + " |")
    if max_rows and len(df) > max_rows:
        lines.append("")
        lines.append("_%d of %d rows shown; the full table is the TSV._"
                     % (max_rows, len(df)))
    return "\n".join(lines) + "\n"


#: why a metric can come back missing, quoted verbatim in the report so the
#: reader does not have to guess whether NA means zero
_MISSING_REASONS = {
    "enc": "fewer than 20 codon-equivalents of synonymous information, which "
           "is Wright's estimator refusing to report a number it cannot support",
    "cai_reference": "no informative codons, or no host reference set",
    "cai_genome": "no host genome codon table",
    "fop": "no host optimal-codon set",
    "tai": "no tRNA gene copy numbers supplied",
    "cufs_host": "no synonymous family observed in both gene and host",
    "rscu_distance_host": "no synonymous family observed in both gene and host",
    "gc3s": "no synonymous codons in the sequence",
    "scuo": "no degenerate families observed",
    "dg_mean": "sequence shorter than the folding window, or ViennaRNA absent",
    "dg_z": "the folding null had fewer than two usable values",
}


#: columns that describe a sequence rather than measure it, so a blank one is
#: missing annotation and not a metric that could not be computed
_NOT_METRICS = ("id", "organism", "taxid", "product", "gene", "source", "note")


def _missing_table(df):
    """Count the metrics that came back missing, with the reason for each.

    A column of nothing but ``None`` arrives as an object column rather than a
    float one, which is exactly the case that matters here: a metric nobody
    could compute for any sequence, such as tAI without tRNA data. Selecting on
    dtype would drop it and the report would claim nothing was missing.
    """
    rows = []
    for col in df.columns:
        if col in _NOT_METRICS:
            continue
        n = int(df[col].isna().sum())
        if n:
            rows.append({"metric": col, "n_missing": n, "n_total": len(df),
                         "reason": _MISSING_REASONS.get(
                             col, "not computed for these sequences")})
    return pd.DataFrame(rows)


def write_report(result, out_dir):
    """Write ``report.md`` and the figures, from the dict :func:`analyse` returns."""
    from pathlib import Path
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    m = result["manifest"]
    per_gene = result["per_gene"]
    host = result.get("host") or {}

    host_metrics = _host_metric_frame(host)
    if len(host_metrics):
        from .pipeline import write_tsv
        write_tsv(host_metrics, out_dir / "host_gene_metrics.tsv")

    figures = []
    try:
        figures.append(("ENC against GC3s with the Wright expected curve",
                        plot_enc_gc3s(per_gene, host_metrics, out_dir)))
        figures.append(("Neutrality plot, GC12 against GC3",
                        plot_neutrality(per_gene, out_dir)))
        figures.append(("PR2 bias at fourfold-degenerate sites",
                        plot_pr2(per_gene, out_dir)))
        rscu_path = out_dir / "rscu_per_gene.tsv"
        if rscu_path.exists():
            rt = pd.read_csv(rscu_path, sep="\t")
            figures.append(("RSCU heatmap, genes against the host",
                            plot_rscu_heatmap(
                                rt, (host.get("table") or {}).get("rscu"),
                                out_dir)))
        figures.append(("Gene against host compatibility",
                        plot_host_compatibility(per_gene, host_metrics, out_dir)))
    except Exception as exc:                                # pragma: no cover
        result.setdefault("notes", []).append("figure generation failed: %s" % exc)

    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    qc_summary = result["qc_summary"]
    n_req, n_ok = m["n_requested"], m["n_passed_qc"]

    L = []
    a = L.append
    a("# codonamr report\n")
    a("Generated %s by codonamr %s. Seed %s, so this run is reproducible with "
      "the same inputs.\n" % (now, m["codonamr_version"], m["seed"]))

    a("\n## What was analysed\n")
    a(_md_table(pd.DataFrame([{
        "sequences requested": n_req,
        "passed QC": n_ok,
        "rejected": n_req - n_ok,
        "host": m["host"] or "none",
        "host assembly": m["host_assembly"] or "none",
        "host CDS": m["host_n_cds"],
        "reference set genes": m["reference_set_size"],
        "permutations per gene": m["n_permutations"],
        "folding": "yes" if m["folding"] else "no",
        "runtime (s)": m["runtime_seconds"],
    }]).T.reset_index().rename(columns={"index": "field", 0: "value"})))

    a("\n## Quality control\n")
    a("Every rejection reason and its count. Nothing was dropped silently.\n")
    a(_md_table(qc_summary))
    rejected = result["qc"][result["qc"]["reason"] != "ok"]
    if len(rejected):
        a("\nRejected sequences:\n")
        a(_md_table(rejected[["id", "reason", "n_nt_input", "organism"]],
                    max_rows=50))

    a("\n## Missing values\n")
    miss = _missing_table(per_gene)
    if len(miss):
        a("A metric reported NA was not computed, and the reason is given. "
          "None of these sequences were removed from the other tables.\n")
        a(_md_table(miss))
    else:
        a("Every metric returned a value for every sequence that passed QC.\n")

    a("\n## Codon usage\n")
    cols = [c for c in ("n_codons", "gc", "gc3s", "enc", "scuo",
                        "cai_reference", "fop", "tai", "cufs_host")
            if c in per_gene]
    desc = per_gene[cols].describe().T.reset_index().rename(
        columns={"index": "metric"})
    a(_md_table(desc[["metric", "count", "mean", "std", "min", "50%", "max"]]))
    a("\nPer-gene values are in `metrics.tsv` and `per_gene_all.tsv`.\n")

    a("\n## Controls\n")
    a("### Length-matched null\n")
    cl = result.get("controls_length_matched")
    if cl is not None and len(cl):
        a("Each gene is compared against host genes of its own length, because "
          "ENC is unstable below roughly 200 codons and a short gene looks "
          "unbiased for reasons that have nothing to do with selection. "
          "`z` is against the length-matched host distribution, not the genome "
          "mean.\n")
        summ = (cl.groupby("metric")
                .agg(n_genes=("id", "nunique"), mean_z=("z", "mean"),
                     n_p_below_0_05=("p", lambda s: int((s < 0.05).sum())),
                     median_pool=("pool_size", "median"),
                     drawn_with_replacement=("with_replacement", "sum"))
                .reset_index())
        a(_md_table(summ))
        if bool(cl["with_replacement"].any()):
            a("\nSome nulls were drawn with replacement because too few host "
              "genes matched on length. Those p-values are lower bounds.\n")
    else:
        a("No length-matched null was built; see the notes below.\n")

    a("\n### Composition control\n")
    cc = result.get("controls_composition")
    if cc is not None and len(cc):
        a("Partial Spearman correlation controlling for the A, C and G "
          "fractions separately (T is dropped because the four sum to 1). "
          "`r_raw` is the uncontrolled correlation. Where the two disagree, "
          "the uncontrolled one was measuring composition.\n")
        a(_md_table(cc))
    else:
        a("Not enough genes with complete values to run a composition "
          "control.\n")

    a("\n### Do the controls work\n")
    val = result.get("validation")
    if val is not None and len(val):
        a("`injection_recovery` plants a partial correlation of known size in "
          "data with a composition confounder and reports what each control "
          "returns. `false_positive_rate` plants nothing and reports how often "
          "an association is found anyway; for a control that works, "
          "`rejection_rate` should be close to 0.05 there. These numbers come "
          "from this dataset's own base composition where there were enough "
          "genes to resample it.\n")
        a(_md_table(val))
    else:
        a("Validation was not run.\n")

    a("\n## Randomisation\n")
    rand = result.get("randomisation")
    if rand is not None and len(rand):
        a("Within-family permutation holds codon composition, ENC, CAI and "
          "GC3s exactly fixed, so anything that moves is attributable to codon "
          "order alone. `dinuc_residual` is the leftover L1 distance in "
          "junction dinucleotide composition after the dinucleotide-controlled "
          "shuffle: 0 means the control succeeded, and a non-zero value means "
          "part of any folding effect could still be dinucleotide composition.\n")
        keep = [c for c in ("junction_chi2_z", "junction_chi2_p",
                            "dinuc_residual_mean", "dinuc_residual_zero_fraction",
                            "dg_mean", "dg_z", "dg_p", "dg_dinuc_z",
                            "dg_dinuc_p", "enc_if_max_bias", "enc_if_uniform")
                if c in rand]
        if keep:
            d = rand[keep].describe().T.reset_index().rename(
                columns={"index": "statistic"})
            a(_md_table(d[["statistic", "count", "mean", "std", "min", "max"]]))
        a("\nPer-gene values are in `randomisation.tsv`.\n")

    a("\n## Figures\n")
    for label, paths in figures:
        if not paths:
            continue
        png, pdf = paths
        a("\n### %s\n" % label)
        a("![%s](figures/%s)\n" % (label, png.name))
        a("\nVector version: `figures/%s`.\n" % pdf.name)

    notes = result.get("notes") or []
    a("\n## Notes and caveats for this run\n")
    if notes:
        for n in notes:
            a("- %s" % n)
    else:
        a("- nothing unusual was recorded during this run")
    a("")
    a("\n## Standing limitations\n")
    a("- Accession counts reflect what was deposited, not what is prevalent in "
      "any population. Nothing here is an epidemiological estimate.")
    a("- CAI is defined against a reference set from one genome. CAI values "
      "computed against different hosts are not on the same scale and should "
      "not be compared directly.")
    a("- ENC on genes under about 200 codons is biased upward and noisy, which "
      "is what the length-matched null is for. Prefer SCUO for short genes.")
    if m["folding"]:
        a("- Folding energies come from a nearest-neighbour model at %s, on "
          "naked RNA, with no ribosomes, no proteins and no ions beyond the "
          "model's assumptions. Treat them as a comparative statistic between "
          "sequences, not as a physical energy."
          % ("%s C" % m["temperature_c"] if m["temperature_c"] is not None
             else "the ViennaRNA default of 37 C"))
    a("")
    (out_dir / "report.md").write_text("\n".join(L))
    return out_dir / "report.md"


def _host_metric_frame(host, max_genes=4000):
    """ENC, GC3s and CAI for the host's own genes, the background distribution.

    Capped at ``max_genes`` because the only use is a background cloud on the
    ENC plot and a histogram, and folding a whole genome's worth of indices for
    that is not a good use of the user's time. The cap is stated in the report.
    """
    seqs = (host or {}).get("cds_clean") or []
    if not seqs:
        return pd.DataFrame(columns=["gc3s", "enc", "cai_reference"])
    from .metrics import cai as _cai, enc as _enc, gc3s as _gc3s
    w = host.get("w_reference")
    step = max(1, len(seqs) // max_genes)
    rows = []
    for s in seqs[::step]:
        rows.append({"gc3s": _gc3s(s), "enc": _enc(s),
                     "cai_reference": _cai(s, w) if w else None})
    return pd.DataFrame(rows)
