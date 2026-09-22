# API reference

The public API of `codonamr`. Everything listed here is importable from the top
level, for example `from codonamr import analyse, enc, cufs`. Formulae and
primary references for the metrics are in [metrics.md](metrics.md).

Conventions used throughout:

- Sequences are plain uppercase DNA strings, in frame, with the terminal stop
  already removed by `check_cds`.
- A metric that cannot be computed returns `None`. It never returns a
  placeholder number, and the pipeline reports the `None` rather than dropping
  the sequence.
- Randomisation functions take a `random.Random` instance so a caller controls
  reproducibility.

---

## Pipeline

### `analyse(ids=None, fasta=None, out_dir="codonamr_results", email=None, api_key=None, host=None, folding=True, n_permutations=100, seed=0, threads=None, ...)`

Runs the whole analysis and writes everything to `out_dir`. Either `ids` or
`fasta` is required. `ids` may be a list of accessions or the path of a file
with one accession per line. `host` overrides the organism inferred from the
records and may be a taxid, an organism name, or the path of a CDS FASTA.
`threads=1` runs everything in the calling process. `seed` fully determines the
output.

Other arguments: `cache_dir`, `db`, `min_len`, `require_stop`, `require_start`,
`trnascan`, `temperature`, `fold_window`, `fold_step`, `n_folding`,
`n_length_matched`, `host_max_cds`, `validate`, `verbose`.

Returns a dict with `manifest`, `qc`, `qc_summary`, `metrics`, `per_gene`,
`randomisation`, `controls_length_matched`, `controls_composition`,
`validation`, `host`, `records`, `notes` and `out_dir`. The DataFrames are the
same tables written to disk.

### `codonamr.report.write_report(result, out_dir)`

Writes `report.md` and the figures from what `analyse` returns. Called by
`analyse`, and usable on its own to regenerate a report.

---

## Sequence input

### `fetch_cds(ids, email, api_key=None, cache_dir=None, db="nuccore", batch_size=200, verbose=True)`

Fetches coding sequences from NCBI. `email` is required by NCBI policy and a
`ValueError` is raised without one. `db="auto"` guesses per accession. Protein
accessions are resolved to their coding nucleotide region through the
`/coded_by` qualifier. Requests are batched 200 at a time, throttled to 3 per
second without an API key and 10 with one, retried with exponential backoff on
HTTP 429 and 5xx, and cached on disk under the accession.

Returns a list of `CdsRecord`. One accession can yield several records if it
carries several CDS features, so do not assume the output has the same length as
the input.

### `read_fasta(path, organism=None, taxid=None)`

Reads local sequences, so the package is usable with no network at all.
`[organism=...]`, `[gene=...]` and `[protein=...]` tags in the description are
parsed if present, which is the format of NCBI CDS downloads.

### `resolve_host_genome(taxid_or_organism, email=None, api_key=None, cache_dir=None, divisions=("bacteria", "archaea"), max_cds=None, verbose=True)`

Finds the RefSeq reference or representative genome for an organism and
downloads its CDS set. This is what makes the host comparison automatic.
Matching falls back from exact taxid to species taxid to exact organism name to
genus, and the match that was used is recorded in the returned `AssemblyHit`.
Returns `{"assembly": AssemblyHit, "cds": [CdsRecord, ...]}`, or `None` when no
reference or representative genome exists for the organism.

### `CdsRecord`

Dataclass with `id`, `sequence`, `description`, `organism`, `taxid`, `product`,
`gene`, `source` and `extra`. `as_dict()` returns it as a plain dict.

---

## Quality control

### `check_cds(seq, min_len=90, require_stop=True, require_start=True, starts=("ATG", "GTG", "TTG"))`

Returns `(cleaned_sequence_or_None, reason)`. The cleaned sequence is uppercase
DNA with the terminal stop removed. `reason` is one of `REASONS`: `ok`,
`too_short`, `not_multiple_of_three`, `ambiguous_bases`, `no_stop_codon`,
`internal_stop`, `bad_start`.

### `clean_cds(seq, **kw)`

Returns the cleaned sequence only, or `None`.

---

## Metrics

| Function | Returns |
|---|---|
| `codon_counts(seq, drop_stops=True)` | counts of every sense codon, including zeros |
| `aa_composition(seq)` | amino-acid counts, stops excluded |
| `rscu(seq_or_counts)` | RSCU per codon, `None` for an unused family |
| `enc(seq, min_expected=20)` | effective number of codons, or `None` if uninformative |
| `build_w(ref_seqs)` | relative adaptiveness `w` from a reference set |
| `cai(seq, w)` | codon adaptation index, or `None` if no informative codon |
| `gc_content(seq)` | overall GC fraction |
| `gc3s(seq)` | GC at synonymous third positions |
| `fop(seq, optimal_codons)` | frequency of optimal codons |
| `scuo(seq)` | synonymous codon usage order, 0 to 1 |

---

## Host comparison

### `host_codon_table(cds_seqs)`

Pooled `counts`, `rscu`, `w`, `n_genes` and `n_codons` for a genome. The `w`
here is genome-wide and is not the CAI `w`, which is defined against a highly
expressed reference set. It is an explicit fallback for genomes where no
reference set can be identified, and the report says when that fallback was
used.

### `reference_set(cds_seqs, products, genes=None, min_genes=20, include_elongation_factors=False, extra_patterns=())`

Picks a highly expressed reference set for CAI from ribosomal protein product
strings. Returns an empty list if fewer than `min_genes` are found, so the
caller can fall back rather than compute a CAI from five genes.

Ribosomal protein *modification enzymes* carry the words "ribosomal protein" in
their product names: prmA and prmB are L11 and L3 methyltransferases, rimK
ligates glutamate onto S6, rimO is an S12 methylthiotransferase, roxA
hydroxylates L16, ycaO is its paralogue. They are ordinary low-abundance
enzymes, and letting them in pulls `w` towards average genome usage in exactly
the families where real ribosomal proteins are most biased. They are excluded by
gene symbol and by product phrase.

### `optimal_codons(rscu_table, threshold=1.0)`

The codons the host over-uses, that is RSCU > 1, as a frozenset.

### `tai(seq, trna_gene_copy_numbers, prokaryote=True, include_singletons=False, weights=None)`

tRNA adaptation index. `trna_gene_copy_numbers` maps anticodon, uppercase DNA,
to gene copy number, with the bacterial lysidine tRNA-Ile under the key
`CAT_Ile2`. Returns `None` when no tRNA data are available: a missing tAI is
reported, never imputed. Pass `weights` from `tai_weights` to avoid rebuilding
the table for every gene.

### `trna_copy_numbers_from_tRNAscan(path, include_pseudo=False, include_undetermined=False)`

Counts tRNA genes per anticodon from tRNAscan-SE output, versions 1.x and 2.x.
Columns are located by shape rather than by fixed offset because the layout
changed between versions.

### `cufs(rscu_a, rscu_b)`

Codon usage frequency similarity, the square root of the Jensen-Shannon
divergence in bits. A distance in [0, 1] despite the name: 0 is identical usage.
Families empty in either table are skipped, and `None` is returned if no family
is shared.

### `rscu_distance(rscu_a, rscu_b, metric="euclidean")`

Distance between two RSCU tables over the families both define. `metric` is one
of `euclidean`, `manhattan`, `cosine`, `correlation`.

---

## Randomisation

### `randomise_replace(seq, preferred, p, rng)`

Replaces each codon with a synonymous one, the family's preferred codon drawn
with marginal probability `p`. Changes codon composition, holds the peptide
fixed.

### `randomise_shuffle(seq, rng)`

Permutes codon positions within each synonymous family. Peptide, codon counts,
ENC, CAI and GC3s are all invariant. Only codon order changes.

### `randomise_shuffle_dinuc(seq, rng, max_iter=None, tol=0)`

Design B with the junction dinucleotide composition held as close to the native
profile as a constrained search can get. Returns `(sequence, residual_L1)`. A
residual of 0 means the match is exact. Report the residual, do not assume it
reached zero.

### `junction_counts(codons)`

Dinucleotides spanning codon boundaries, base 3 of codon `i` with base 1 of
codon `i+1`. These are the only part of dinucleotide composition that codon
order can change.

---

## Controls

### `length_matched_null(seq, host_cds, n=100, rng=None, tolerance=0.2, min_pool=20, widen=(0.2, 0.35, 0.5, 1.0))`

Draws `n` host genes of the query gene's codon length. Returns a
`LengthMatchedNull` with `sequences`, `target_codons`, `pool_size`, `tolerance`,
`with_replacement` and `note`. Iterating over it yields the sequences.

ENC is unstable below roughly 200 codons, so a short gene looks unbiased for
reasons unrelated to selection. Matching on length makes the bias apply to both
sides of the comparison. If the pool is too small even after widening, genes are
drawn with replacement and `with_replacement` records it, which makes any
p-value from that null a lower bound.

### `empirical_z(value, null_values)`

Returns `z`, `p`, `n_null`, `null_mean` and `null_sd`. The p-value uses the
`(r + 1) / (m + 1)` correction, so it is never exactly zero.

### `base_fractions(seq)`

The fractions of A, C, G and T as a 4-tuple.

### `spearman_partial(x, y, covariates)`

Spearman partial correlation. Returns `r`, `p`, `n`, `k` and `df`. Removes
monotone dependence on the covariates, not arbitrary dependence.

### `composition_control(values, gc_per_nucleotide, response, method="spearman", drop_base=3)`

Correlates `values` with `response`, controlling for base composition.
`gc_per_nucleotide` is an `(n, 4)` array of A, C, G, T fractions; a 1-D array is
accepted and reproduces the usual GC-only control for comparison. `drop_base=3`
drops the T column because the four fractions sum to 1 and the design matrix
would otherwise be singular. Returns `r_raw`, `p_raw`, `r_partial`, `p_partial`,
`n`, `k` and `df`.

### `injection_recovery(n=200, effect=0.3, confound=0.8, n_replicates=200, seed=0, alpha=0.05, sequences=None, skew=0.5)`

Plants a partial correlation of exactly `effect` between two variables that also
share a composition confounder, then analyses the data with no control, with a
GC-only control, and with the per nucleotide composition control. `skew` is the
fraction of the confounder carried by AT skew, which a GC-only control cannot
see. `sequences` is an optional `(m, 4)` array of real base fractions to
resample, which makes the validation specific to a dataset.

Returns a DataFrame with one row per control: `mean_r`, `sd_r`, `bias` and
`rejection_rate`.

### `false_positive_rate(n=200, confound=0.8, n_replicates=500, seed=0, alpha=0.05, sequences=None, skew=0.5)`

The same design with `effect=0`, so `rejection_rate` is the false positive rate.
It should be close to `alpha` for a control that works.

---

## Folding

### `folding_available()`

True if ViennaRNA can be imported.

### `dg_windowed(seq, window=40, step=1, temperature=None)`

Mean local folding free energy over sliding windows. `temperature` in Celsius
sets the ViennaRNA model temperature. Returns `None` if the sequence is shorter
than the window. Raises `ImportError` with an actionable message if ViennaRNA is
absent.

---

## Figures

All in `codonamr.report`, all returning `(png_path, pdf_path)` when given an
`out_dir` and a matplotlib figure otherwise: `plot_enc_gc3s`,
`plot_neutrality`, `plot_pr2`, `plot_rscu_heatmap`, `plot_host_compatibility`.
`wright_expected_enc(gc3s)` gives the Wright expected ENC curve. `PALETTE` is
the Okabe and Ito colourblind-safe palette used throughout.

---

## Genetic code

`CODONS` in TCAG order, `CODON2AA`, `FAMILY` (amino acid to its codons),
`DEGENERATE` (families with more than one codon), and `codon_list(seq)`, which
splits an in-frame DNA string into codons.
