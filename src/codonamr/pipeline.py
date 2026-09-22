"""End to end analysis: identifiers in, tables and a report out.

The design goal is that a user supplies sequence identifiers and nothing else.
Everything the comparison needs, including which genome to compare against, is
derived from the records themselves. Every intermediate is written as a TSV, so
no result is only available from inside the program, and every stochastic step
is seeded from the record identifier rather than from call order, which is what
makes the output identical whether it was produced on one process or on twelve.

Order of operations, which is also the order of the output files: fetch, QC,
host resolution, per-gene metrics, host comparison, controls, randomisation,
folding, report.
"""
import hashlib
import json
import multiprocessing as mp
import os
import random
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd

from . import folding as folding_mod
from . import qc
from .controls import (base_fractions, composition_control, empirical_z,
                       false_positive_rate, injection_recovery,
                       length_matched_null)
from .fetch import (CdsRecord, default_cache_dir, fetch_cds, read_fasta,
                    resolve_host_genome)
from .genetic_code import CODON2AA, CODONS, FAMILY, codon_list
from .host import (cufs, host_codon_table, optimal_codons, reference_set,
                   rscu_distance, tai, tai_weights,
                   trna_copy_numbers_from_tRNAscan)
from .metrics import (aa_composition, cai, codon_counts, enc, fop, gc3s,
                      gc_content, rscu, scuo)
from .randomise import (junction_counts, randomise_replace, randomise_shuffle,
                        randomise_shuffle_dinuc)

__all__ = ["analyse", "write_tsv"]

VERSION = "0.1.0"


# ----------------------------------------------------------------- helpers
def _seed_for(seed, key):
    """Per-record seed, stable across runs, processes and input order.

    ``hash()`` is salted per process in Python 3, so it cannot be used here.
    A digest of the record id keeps every gene's randomisation reproducible
    even if the gene list is reordered or the work is split differently.
    """
    h = hashlib.sha256(("%s|%s" % (seed, key)).encode()).digest()
    return int.from_bytes(h[:8], "big")


def write_tsv(df, path, index=False):
    """Write a table, creating parent directories. Every table goes through here."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index, float_format="%.6g", na_rep="NA")
    return path


def _map(func, items, threads):
    """Parallel map that keeps input order, with an honest single-process path.

    ``threads=1`` never touches multiprocessing at all, so the package stays
    debuggable and usable inside environments that forbid process creation.
    """
    items = list(items)
    if not items:
        return []
    if threads is None:
        threads = min(8, os.cpu_count() or 1)
    threads = max(1, int(threads))
    if threads == 1 or len(items) == 1:
        return [func(x) for x in items]
    try:
        ctx = mp.get_context("fork" if hasattr(os, "fork") else "spawn")
        with ctx.Pool(threads) as pool:
            return pool.map(func, items, chunksize=1)
    except (OSError, ValueError, RuntimeError) as exc:
        print("codonamr: multiprocessing unavailable (%s), running serially"
              % exc, file=sys.stderr)
        return [func(x) for x in items]


# ------------------------------------------------------------- input stage
def _load_records(ids, fasta, email, api_key, cache_dir, db, verbose):
    records = []
    if fasta:
        for p in ([fasta] if isinstance(fasta, (str, Path)) else list(fasta)):
            records.extend(read_fasta(p))
    if ids:
        wanted = list(ids)
        if len(wanted) == 1 and Path(str(wanted[0])).exists():
            wanted = Path(str(wanted[0])).read_text().split()
        records.extend(fetch_cds(wanted, email=email, api_key=api_key,
                                 cache_dir=cache_dir, db=db, verbose=verbose))
    seen, out = set(), []
    for r in records:
        if r.id in seen:
            continue
        seen.add(r.id)
        out.append(r)
    return out


def _qc_stage(records, min_len, require_stop, require_start):
    rows, kept = [], []
    for r in records:
        clean, reason = qc.check_cds(r.sequence, min_len=min_len,
                                     require_stop=require_stop,
                                     require_start=require_start)
        rows.append({"id": r.id, "reason": reason, "n_nt_input": len(r.sequence),
                     "n_codons_kept": len(codon_list(clean)) if clean else 0,
                     "organism": r.organism, "taxid": r.taxid,
                     "product": r.product, "gene": r.gene, "source": r.source})
        if clean:
            kept.append(CdsRecord(id=r.id, sequence=clean,
                                  description=r.description,
                                  organism=r.organism, taxid=r.taxid,
                                  product=r.product, gene=r.gene,
                                  source=r.source, extra=dict(r.extra)))
    return kept, pd.DataFrame(rows)


# ----------------------------------------------------------- metric stage
def _metric_job(job):
    """One gene's metrics. Module level and dict-argumented so it can be pickled."""
    seq = job["sequence"]
    w_ref = job.get("w_reference")
    w_gen = job.get("w_genome")
    opt = set(job.get("optimal_codons") or ())
    host_rscu = job.get("host_rscu")
    tai_w = job.get("tai_weights")
    r = rscu(seq)
    a, c, g, t = base_fractions(seq)
    out = {
        "id": job["id"],
        "n_codons": len(codon_list(seq)),
        "gc": gc_content(seq),
        "gc3s": gc3s(seq),
        "frac_A": a, "frac_C": c, "frac_G": g, "frac_T": t,
        "enc": enc(seq),
        "scuo": scuo(seq),
        "cai_reference": cai(seq, w_ref) if w_ref else None,
        "cai_genome": cai(seq, w_gen) if w_gen else None,
        "fop": fop(seq, opt) if opt else None,
        "tai": tai(seq, None, weights=tai_w) if tai_w else None,
        "cufs_host": cufs(r, host_rscu) if host_rscu else None,
        "rscu_distance_host": rscu_distance(r, host_rscu) if host_rscu else None,
    }
    out["_rscu"] = r
    out["_aa"] = aa_composition(seq)
    out["_counts"] = codon_counts(seq)
    return out


def _gc12_gc3(seq):
    """Neutrality-plot coordinates: mean GC at positions 1 and 2, and GC3.

    All codons are used here, including Met, Trp and the third positions of
    single-codon families, because the neutrality plot is a statement about
    mutational pressure across the gene rather than about synonymous choice.
    This is why GC3 here differs slightly from the GC3s used for the ENC plot.
    """
    g1 = g2 = g3 = n = 0
    for c in codon_list(seq):
        if CODON2AA.get(c) is None:
            continue
        n += 1
        g1 += c[0] in "GC"
        g2 += c[1] in "GC"
        g3 += c[2] in "GC"
    if not n:
        return None, None
    return (g1 + g2) / (2.0 * n), g3 / float(n)


def _pr2(seq):
    """PR2 coordinates at fourfold-degenerate sites: A3/(A3+T3) and G3/(G3+C3)."""
    a = t = g = c = 0
    for cod in codon_list(seq):
        aa = CODON2AA.get(cod)
        if aa is None or aa == "*" or len(FAMILY[aa]) != 4:
            continue
        b = cod[2]
        a += b == "A"
        t += b == "T"
        g += b == "G"
        c += b == "C"
    at = (a / (a + t)) if (a + t) else None
    gc = (g / (g + c)) if (g + c) else None
    return at, gc


# ------------------------------------------------------- randomisation job
def _randomise_job(job):
    """Permutation test on codon order for one gene.

    Design B (within-family permutation) holds codon composition, and therefore
    ENC, CAI and GC3s, exactly fixed, so anything it moves is attributable to
    codon order. Two statistics are tested. The first is the deviation of
    junction dinucleotide counts from their permutation expectation, which needs
    no external software. The second is mean local folding energy, which needs
    ViennaRNA and is only computed when it is available and asked for.
    """
    seq = job["sequence"]
    rng = random.Random(job["seed"])
    n = int(job["n_permutations"])
    out = {"id": job["id"], "n_permutations": n}

    obs_j = junction_counts(codon_list(seq))
    keys = sorted(set(obs_j) | set("".join(p) for p in
                                   [(x, y) for x in "ACGT" for y in "ACGT"]))
    null_j = []
    shuffled = []
    for _ in range(n):
        s = randomise_shuffle(seq, rng)
        shuffled.append(s)
        null_j.append(junction_counts(codon_list(s)))
    if null_j:
        mean = {k: sum(d.get(k, 0) for d in null_j) / float(n) for k in keys}

        def chi2(d):
            return sum((d.get(k, 0) - mean[k]) ** 2 / mean[k]
                       for k in keys if mean[k] > 0)

        res = empirical_z(chi2(obs_j), [chi2(d) for d in null_j])
        out.update({"junction_chi2": chi2(obs_j),
                    "junction_chi2_null_mean": res["null_mean"],
                    "junction_chi2_z": res["z"], "junction_chi2_p": res["p"]})

    n_dinuc = min(n, int(job.get("n_dinuc", 20)))
    residuals = []
    dinuc_seqs = []
    for _ in range(n_dinuc):
        s, resid = randomise_shuffle_dinuc(seq, rng)
        residuals.append(resid)
        dinuc_seqs.append(s)
    if residuals:
        out["dinuc_residual_mean"] = sum(residuals) / len(residuals)
        out["dinuc_residual_max"] = max(residuals)
        out["dinuc_residual_zero_fraction"] = sum(r == 0 for r in residuals) / len(residuals)

    if job.get("folding"):
        window, step = job["fold_window"], job["fold_step"]
        temp = job.get("temperature")
        obs = folding_mod.dg_windowed(seq, window=window, step=step,
                                      temperature=temp)
        out["dg_mean"] = obs
        k = min(len(shuffled), int(job.get("n_folding", 20)))
        null = [folding_mod.dg_windowed(s, window=window, step=step,
                                        temperature=temp)
                for s in shuffled[:k]]
        res = empirical_z(obs, null)
        out.update({"dg_null_mean": res["null_mean"], "dg_null_sd": res["null_sd"],
                    "dg_z": res["z"], "dg_p": res["p"], "dg_n_null": res["n_null"]})
        kd = min(len(dinuc_seqs), k)
        if kd:
            nulld = [folding_mod.dg_windowed(s, window=window, step=step,
                                             temperature=temp)
                     for s in dinuc_seqs[:kd]]
            resd = empirical_z(obs, nulld)
            out.update({"dg_dinuc_null_mean": resd["null_mean"],
                        "dg_dinuc_z": resd["z"], "dg_dinuc_p": resd["p"],
                        "dg_dinuc_n_null": resd["n_null"]})

    pref = job.get("preferred")
    if pref:
        maxbias = randomise_replace(seq, pref, 1.0, rng)
        uniform = randomise_replace(seq, pref, 0.0, rng)
        out["enc_if_max_bias"] = enc(maxbias)
        out["enc_if_uniform"] = enc(uniform)
        out["gc3s_if_max_bias"] = gc3s(maxbias)
        out["gc3s_if_uniform"] = gc3s(uniform)
    return out


# --------------------------------------------------------- control stage
def _control_job(job):
    """Length-matched null for one gene, on the indices that need one."""
    seq = job["sequence"]
    rng = random.Random(job["seed"])
    matched = length_matched_null(seq, job["host_cds"], n=job["n"], rng=rng)
    rows = []
    funcs = {"enc": enc, "scuo": scuo, "gc3s": gc3s}
    w = job.get("w_reference")
    if w:
        funcs["cai_reference"] = lambda s, _w=w: cai(s, _w)
    for name, f in funcs.items():
        obs = f(seq)
        null = [f(s) for s in matched.sequences]
        res = empirical_z(obs, null)
        rows.append({"id": job["id"], "metric": name, "observed": obs,
                     "null_mean": res["null_mean"], "null_sd": res["null_sd"],
                     "z": res["z"], "p": res["p"], "n_null": res["n_null"],
                     "n_codons": matched.target_codons,
                     "pool_size": matched.pool_size,
                     "length_tolerance": matched.tolerance,
                     "with_replacement": matched.with_replacement,
                     "note": matched.note})
    return rows


# ----------------------------------------------------------------- driver
def analyse(ids=None, fasta=None, out_dir="codonamr_results", email=None,
            api_key=None, host=None, folding=True, n_permutations=100,
            seed=0, threads=None, cache_dir=None, db="auto",
            min_len=90, require_stop=True, require_start=False,
            trnascan=None, temperature=None, fold_window=40, fold_step=3,
            n_folding=20, n_length_matched=100, host_max_cds=None,
            validate=True, verbose=True):
    """Run the whole analysis and write everything to ``out_dir``.

    Parameters
    ----------
    ids : iterable of str or path, optional
        NCBI accessions, or the path of a file with one accession per line.
    fasta : path or iterable of path, optional
        Local FASTA files, used instead of or alongside ``ids``. With ``fasta``
        alone the run makes no network calls except host resolution, and none at
        all if ``host`` is also supplied as a local file.
    out_dir : path
        Output directory. Created if absent.
    email : str, optional
        Required by NCBI for any network access. Falls back to the
        ``CODONAMR_EMAIL`` environment variable.
    api_key : str, optional
        NCBI API key, or the ``CODONAMR_API_KEY`` environment variable.
    host : str or path, optional
        Taxid or organism name to use as the host, overriding the organism
        inferred from the records, or the path of a FASTA file of host CDS.
    folding : bool
        Compute local folding energies. Silently downgraded with a recorded
        note if ViennaRNA is not installed.
    n_permutations : int
        Permutations per gene for the codon-order test.
    seed : int
        Master seed. Per-gene seeds are derived from it and the gene id, so the
        result does not depend on ``threads`` or on the order of the input.
    threads : int, optional
        Worker processes. ``1`` runs everything in this process. ``None`` uses
        up to 8 cores.
    require_start : bool
        Off by default, because a fetched CDS region is often correct without
        beginning at an annotated start codon, and rejecting those loses real
        genes. QC still reports how many would have failed.
    validate : bool
        Run :func:`codonamr.controls.injection_recovery` and
        :func:`codonamr.controls.false_positive_rate` on the composition of
        this dataset and write the result. On by default: the controls are part
        of the result, not a separate exercise.

    Returns
    -------
    dict
        Tables, the host record and a run manifest. The same contents are on
        disk as TSV.
    """
    t0 = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    email = email or os.environ.get("CODONAMR_EMAIL")
    api_key = api_key or os.environ.get("CODONAMR_API_KEY")
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    notes = []

    # ---- input -----------------------------------------------------------
    records = _load_records(ids, fasta, email, api_key, cache_dir, db, verbose)
    if not records:
        raise ValueError("no sequences: pass ids=..., fasta=..., or both")

    kept, qc_table = _qc_stage(records, min_len, require_stop, require_start)
    write_tsv(qc_table, out_dir / "qc.tsv")
    qc_summary = (qc_table.groupby("reason").size().rename("n_sequences")
                  .reset_index().sort_values("n_sequences", ascending=False))
    write_tsv(qc_summary, out_dir / "qc_summary.tsv")
    if verbose:
        print("codonamr: %d/%d sequences passed QC"
              % (len(kept), len(records)), file=sys.stderr)
    if not kept:
        raise ValueError("every sequence failed QC, see qc_summary.tsv")

    # ---- host ------------------------------------------------------------
    host_info = _resolve_host(kept, host, email, api_key, cache_dir,
                              host_max_cds, min_len, notes, verbose)
    host_cds = host_info["cds_clean"]
    host_table = host_info["table"]
    ref_w = host_info["w_reference"]
    opt = host_info["optimal_codons"]

    trna = None
    tai_w = None
    if trnascan:
        trna = trna_copy_numbers_from_tRNAscan(trnascan)
        tai_w = tai_weights(trna)
        if not tai_w:
            notes.append("tRNAscan file %s yielded no usable anticodons, tAI "
                         "was not computed" % trnascan)
    else:
        notes.append("no tRNA gene copy numbers supplied (--trnascan), so tAI "
                     "is reported as NA. tAI needs tRNAscan-SE output for the "
                     "host genome and cannot be inferred from sequence alone.")

    # ---- per-gene metrics ------------------------------------------------
    jobs = [{"id": r.id, "sequence": r.sequence, "w_reference": ref_w,
             "w_genome": host_table["w"] if host_table else None,
             "optimal_codons": sorted(opt) if opt else None,
             "host_rscu": host_table["rscu"] if host_table else None,
             "tai_weights": tai_w} for r in kept]
    metric_rows = _map(_metric_job, jobs, threads)

    rscu_rows, aa_rows, count_rows = [], [], []
    for row, rec in zip(metric_rows, kept):
        rscu_rows.append(dict({"id": row["id"]}, **row.pop("_rscu")))
        aa_rows.append(dict({"id": row["id"]}, **row.pop("_aa")))
        count_rows.append(dict({"id": row["id"]}, **row.pop("_counts")))
        gc12, gc3 = _gc12_gc3(rec.sequence)
        at3, gc3_pr2 = _pr2(rec.sequence)
        row.update({"gc12": gc12, "gc3": gc3, "pr2_a3_at3": at3,
                    "pr2_g3_gc3": gc3_pr2, "organism": rec.organism,
                    "taxid": rec.taxid, "product": rec.product,
                    "gene": rec.gene, "source": rec.source})
    metrics = pd.DataFrame(metric_rows)
    write_tsv(metrics, out_dir / "metrics.tsv")
    write_tsv(pd.DataFrame(rscu_rows), out_dir / "rscu_per_gene.tsv")
    write_tsv(pd.DataFrame(aa_rows).fillna(0), out_dir / "aa_composition.tsv")
    write_tsv(pd.DataFrame(count_rows), out_dir / "codon_counts.tsv")

    if host_table:
        write_tsv(_host_table_frame(host_table, ref_w, opt),
                  out_dir / "host_codon_table.tsv")

    # ---- controls --------------------------------------------------------
    controls_len = pd.DataFrame()
    if host_cds:
        cjobs = [{"id": r.id, "sequence": r.sequence, "host_cds": host_cds,
                  "n": n_length_matched, "w_reference": ref_w,
                  "seed": _seed_for(seed, "lmn:" + r.id)} for r in kept]
        rows = []
        for res in _map(_control_job, cjobs, threads):
            rows.extend(res)
        controls_len = pd.DataFrame(rows)
        write_tsv(controls_len, out_dir / "controls_length_matched.tsv")
    else:
        notes.append("no host CDS set, so no length-matched null was built. "
                     "ENC and CAI are reported without a length control, which "
                     "for genes under 200 codons is not a safe comparison.")

    # ---- randomisation and folding --------------------------------------
    can_fold = bool(folding) and folding_mod.available()
    if folding and not can_fold:
        notes.append("ViennaRNA is not installed, so folding energies were not "
                     "computed. Install the optional extra with "
                     "`pip install codonamr[folding]`.")
    preferred = host_info.get("preferred_codon")
    rjobs = [{"id": r.id, "sequence": r.sequence,
              "n_permutations": n_permutations, "n_folding": n_folding,
              "folding": can_fold, "fold_window": fold_window,
              "fold_step": fold_step, "temperature": temperature,
              "preferred": preferred,
              "seed": _seed_for(seed, "rand:" + r.id)} for r in kept]
    rand = pd.DataFrame(_map(_randomise_job, rjobs, threads))
    write_tsv(rand, out_dir / "randomisation.tsv")

    merged = metrics.merge(rand, on="id", how="left")
    write_tsv(merged, out_dir / "per_gene_all.tsv")

    # ---- composition control --------------------------------------------
    comp_rows = []
    comp = merged[["frac_A", "frac_C", "frac_G", "frac_T"]].to_numpy(float)
    pairs = [("gc3s", "enc"), ("enc", "cai_reference"), ("enc", "scuo"),
             ("cai_reference", "cufs_host"), ("enc", "dg_mean"),
             ("cai_reference", "dg_mean"), ("gc3s", "dg_mean")]
    for a, b in pairs:
        if a not in merged or b not in merged:
            continue
        va, vb = merged[a], merged[b]
        if va.notna().sum() < 8 or vb.notna().sum() < 8:
            continue
        res = composition_control(va.to_numpy(float), comp, vb.to_numpy(float))
        comp_rows.append(dict({"x": a, "y": b}, **res))
    controls_comp = pd.DataFrame(comp_rows)
    if len(controls_comp):
        write_tsv(controls_comp, out_dir / "controls_composition.tsv")
    else:
        notes.append("too few genes with complete values for a composition "
                     "control; at least 8 are needed and more like 30 before "
                     "a partial correlation means much.")

    # ---- validation of the controls themselves ---------------------------
    validation = pd.DataFrame()
    if validate:
        comp_ok = comp[~pd.isna(comp).any(axis=1)] if len(comp) else comp
        seqs = comp_ok if len(comp_ok) >= 10 else None
        rec = injection_recovery(n=max(len(merged), 50), n_replicates=200,
                                 seed=seed, sequences=seqs)
        rec.insert(0, "experiment", "injection_recovery")
        fpr = false_positive_rate(n=max(len(merged), 50), n_replicates=200,
                                  seed=seed + 1, sequences=seqs)
        fpr.insert(0, "experiment", "false_positive_rate")
        validation = pd.concat([rec, fpr], ignore_index=True)
        write_tsv(validation, out_dir / "validation.tsv")

    manifest = {
        "codonamr_version": VERSION,
        "seed": seed,
        "threads": threads,
        "n_requested": len(records),
        "n_passed_qc": len(kept),
        "n_permutations": n_permutations,
        "folding": can_fold,
        "fold_window": fold_window,
        "fold_step": fold_step,
        "temperature_c": temperature,
        "host": host_info["label"],
        "host_assembly": host_info["accession"],
        "host_match": host_info["match"],
        "host_n_cds": len(host_cds),
        "reference_set_size": host_info["reference_set_size"],
        "tai_available": bool(tai_w),
        "cache_dir": str(cache_dir),
        "runtime_seconds": round(time.time() - t0, 2),
        "notes": notes,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    result = {"manifest": manifest, "qc": qc_table, "qc_summary": qc_summary,
              "metrics": metrics, "per_gene": merged,
              "randomisation": rand, "controls_length_matched": controls_len,
              "controls_composition": controls_comp, "validation": validation,
              "host": host_info, "records": kept, "notes": notes,
              "out_dir": str(out_dir)}

    from .report import write_report
    write_report(result, out_dir)
    return result


def _host_table_frame(host_table, ref_w, opt):
    rows = []
    for c in CODONS:
        aa = CODON2AA[c]
        if aa == "*":
            continue
        rows.append({"codon": c, "aa": aa, "family_size": len(FAMILY[aa]),
                     "count": host_table["counts"].get(c, 0),
                     "rscu": host_table["rscu"].get(c),
                     "w_genome": host_table["w"].get(c),
                     "w_reference": ref_w.get(c) if ref_w else None,
                     "optimal": c in opt if opt else None})
    return pd.DataFrame(rows)


def _resolve_host(kept, host, email, api_key, cache_dir, host_max_cds, min_len,
                  notes, verbose):
    """Decide which genome the genes are compared against, and build its tables.

    The organism is taken from the records themselves, by majority taxid, so
    the user does not have to know it. If the genes come from several organisms
    the majority wins and the rest are compared against a host that is not
    theirs, which is recorded in the manifest and stated in the report. That
    situation is common with resistance genes and is a real limitation, not a
    detail: a comparison against the wrong host is not evidence of anything.
    """
    info = {"label": None, "accession": None, "match": None, "cds": [],
            "cds_clean": [], "table": None, "w_reference": None,
            "optimal_codons": None, "reference_set_size": 0,
            "preferred_codon": None, "products": [], "genes": []}

    host_records = None
    if host and Path(str(host)).exists():
        host_records = read_fasta(host)
        info["label"] = str(host)
        info["match"] = "local FASTA"
    else:
        key = host
        if key is None:
            taxids = Counter(r.taxid for r in kept if r.taxid)
            organisms = Counter(r.organism for r in kept if r.organism)
            key = (taxids.most_common(1)[0][0] if taxids else
                   organisms.most_common(1)[0][0] if organisms else None)
            if taxids and len(taxids) > 1:
                notes.append("genes come from %d different taxa; the host is "
                             "the commonest one (%s, %d/%d genes) and the rest "
                             "are compared against a genome that is not theirs."
                             % (len(taxids), taxids.most_common(1)[0][0],
                                taxids.most_common(1)[0][1], len(kept)))
        if key is None:
            notes.append("no organism or taxid on any record, so no host genome "
                         "could be resolved. CAI, Fop and the length-matched "
                         "null are unavailable. Supply --host.")
            return info
        try:
            got = resolve_host_genome(key, email=email, api_key=api_key,
                                      cache_dir=cache_dir,
                                      max_cds=host_max_cds, verbose=verbose)
        except Exception as exc:
            notes.append("host genome lookup for %s failed: %s" % (key, exc))
            return info
        if not got:
            notes.append("no RefSeq reference or representative genome for %s, "
                         "so no host comparison was made. Supply --host with a "
                         "taxid, an organism name or a local CDS FASTA." % key)
            return info
        host_records = got["cds"]
        asm = got["assembly"]
        info["label"] = asm.organism
        info["accession"] = asm.accession
        info["match"] = "%s, %s" % (asm.category, asm.level)
        if str(asm.taxid) != str(key) and asm.organism.lower() != str(key).lower():
            notes.append("host genome %s (%s) is not an exact taxid match for "
                         "%s; the codon table is from a relative."
                         % (asm.accession, asm.organism, key))

    clean, products, genes = [], [], []
    for r in host_records:
        s, reason = qc.check_cds(r.sequence, min_len=min_len,
                                 require_stop=False, require_start=False)
        if s:
            clean.append(s)
            products.append(r.product)
            genes.append(r.gene)
    info["cds"] = host_records
    info["cds_clean"] = clean
    info["products"] = products
    info["genes"] = genes
    if not clean:
        notes.append("host CDS set had no sequences that passed QC")
        return info

    info["table"] = host_codon_table(clean)
    info["optimal_codons"] = optimal_codons(info["table"]["rscu"])
    ref = reference_set(clean, products, genes)
    info["reference_set_size"] = len(ref)
    if ref:
        from .metrics import build_w
        info["w_reference"] = build_w(ref)
    else:
        notes.append("fewer than 20 ribosomal-protein genes could be identified "
                     "in the host annotation, so CAI was computed against a "
                     "genome-wide w instead of a highly expressed reference "
                     "set. That is a weaker statistic and is not comparable "
                     "with a published CAI.")
        info["w_reference"] = info["table"]["w"]

    pref = {}
    for aa, fam in FAMILY.items():
        best = max(fam, key=lambda c: info["table"]["counts"].get(c, 0))
        pref[aa] = best
    info["preferred_codon"] = pref
    return info
