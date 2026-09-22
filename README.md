# codonamr: codon usage in antimicrobial resistance genes, with controls that are themselves validated

A Python package for codon usage analysis of antimicrobial resistance genes.
You supply sequence identifiers. It fetches the coding sequences, works out
which organism they came from, resolves that organism to a RefSeq reference
genome, builds the host codon table from it, computes the standard codon usage
indices, and then applies the two controls that decide whether any difference it
finds is real: a length-matched null, because the indices are biased on short
genes and resistance genes are short, and a partial correlation on per
nucleotide composition, because base composition drives codon choice on its own.
The controls are validated inside the package by injecting an association of
known size and by injecting none.

This repository contains **code only**. No sequence data, results tables,
figures or manuscript files are included, every input is downloaded from public
databases by the package itself and written to a directory you name (see
[section 9](#9-output-files)).

---

## 1. Disclosure: use of artificial intelligence

**Artificial intelligence was used in producing this software, and this section
states exactly how.**

This package was written with the assistance of **Claude (Anthropic)**, used
through Claude Code as a programming assistant. That assistance covered:

- drafting and refactoring the modules collected in `src/codonamr/`;
- implementing published algorithms in Python after checking the primary
  literature for each one: the codon adaptation index (Sharp and Li 1987), the
  effective number of codons (Wright 1990), relative synonymous codon usage
  (Sharp, Tuohy and Mosurski 1986), the synonymous codon usage order (Wan et al.
  2004), the tRNA adaptation index with the wobble constants of dos Reis, Savva
  and Wernisch (2004), and the dinucleotide-preserving synonymous shuffle of
  Katz and Burge (2003);
- writing the test suite in `tests/` and the documentation in `docs/`.

**What AI did not do.** No sequence data was generated, altered or imputed by
AI. No result, statistic or figure value produced by this package comes from a
language model: every number a run reports is computed by the code in this
repository from records downloaded from public databases at run time. No
references were invented. Every method implemented here is attributed to its
primary publication, and those attributions are collected in
[`docs/metrics.md`](docs/metrics.md) so that a reader can check each one against
the source.

**Accountability.** The author has reviewed the code and retains full
responsibility for its correctness, for the analyses it performs, and for all
conclusions drawn from it. AI assistance does not transfer any part of that
responsibility.

This disclosure is made because reproducibility depends on knowing how research
software was built, and because this package implements published methods
in house rather than wrapping established binaries. Readers are encouraged to
check those implementations, and
[section 7](#7-validating-the-implementations) says how.

---

## 2. Statement of need

Codon usage indices are easy to compute and easy to over-read. Three problems
recur in the literature on mobile resistance genes.

The first is that the indices are small-sample statistics. The effective number
of codons is estimated from the homozygosity of each synonymous family, and that
estimate is biased and noisy when a family holds only a few codons. Resistance
genes are short, often under 300 codons, so a resistance gene will look less
biased than the chromosomal genes around it for reasons that have nothing to do
with selection. Comparing it against a genome mean makes that artefact look like
a result.

The second is confounding by base composition. GC content drives synonymous
codon choice mutationally, and a gene acquired from a compositionally different
genome differs from its new host in more than its GC total. Controlling for GC
as a single number leaves the strand-asymmetric part of composition, the part
that separates A from T and G from C, entirely uncontrolled.

The third is that the host side of the comparison is usually assembled by hand.
Deciding which genome to compare against, downloading its coding sequences,
identifying the highly expressed genes that define a CAI reference set and
keeping ribosomal protein modification enzymes out of that set are all manual
steps, and manual steps are where analyses diverge between papers.

**Who this is for.** Researchers studying horizontally transferred genes,
resistance gene expression and host adaptation, who want the standard indices
computed the same way every time, with the host comparison built automatically
from the record's own taxonomy, and with the confounders addressed explicitly
rather than acknowledged in a discussion section.

**How it differs from existing tools.** The established tools are good at what
they do and this package does not replace them.

| Tool | What it does well | What it does not do |
|---|---|---|
| **CodonW** (Peden 1999) | The reference implementation of ENC, CAI, Fop, CBI and correspondence analysis. Fast, stable, and the source of the conventions this package follows. | It is a calculator. It does not fetch sequences, resolve a host genome, build a reference set, or provide any null model or statistical control. |
| **CodonW-derived web servers** | Convenient for a handful of sequences with no installation. | Not scriptable or reproducible in a pipeline, usually capped in input size, and the version actually running is often unstated. |
| **CAIcal** (Puigbo et al. 2008) | CAI with an expected-CAI null that corrects for base composition, which is a genuine statistical control and the closest existing work to this package. | The control applies to CAI only. There is no length matching, no per nucleotide composition control, no host genome resolution and no validation of the control itself. |
| **EMBOSS cusp and chips** | Solid, well-tested codon usage tables and ENC within a large and widely installed suite. | Per-file calculators with no host comparison, no controls and no reference set logic. |
| **codonW-style R and Python ports** (for example coRdon, CAI, python-codon-tables) | Convenient library access to the indices inside an existing analysis. | Same scope as the calculators: indices, not designs. |

What is new here is not an index. It is the surrounding design: identifier in,
host resolved automatically, reference set built with the ribosomal protein
modification enzymes deliberately excluded, every index reported with a
length-matched null and a per nucleotide composition control, and both controls
checked by injection recovery and false positive rate measurement that ship as
part of the package and can be run in one command.

---

## 3. Install

From a clone of this repository, which is the way to install it until the
first release is on PyPI:

```bash
git clone https://github.com/bioforkgithub/codonamr
cd codonamr
pip install -e ".[dev]"
```

Once released:

```bash
pip install codonamr
```

Python 3.9 or newer. The required dependencies are biopython, numpy, scipy,
pandas and matplotlib.

Local mRNA folding is optional because ViennaRNA is a larger install than the
rest put together. Everything else works without it, and a run without it says
so in the report rather than silently omitting the folding section:

```bash
pip install "codonamr[folding]"
```

### Configuration

NCBI asks that every E-utilities request carries a contact address, so set one
once and the commands below need no further arguments:

```bash
export CODONAMR_EMAIL="you@example.org"
```

Two optional variables:

| variable | what it does | default |
|---|---|---|
| `CODONAMR_EMAIL` | contact address sent with every NCBI request | none, and `--ids` will refuse to run without it |
| `CODONAMR_API_KEY` | NCBI API key, which raises the rate limit from 3 to 10 requests per second | none |
| `CODONAMR_CACHE` | where downloads are cached | the per-user cache directory |

Downloads are cached on disk under the accession, so a second run of the same
analysis makes no network calls at all. Delete the cache directory to force a
refresh.

---

## 4. Quickstart

Identifiers in, everything out. This is the whole interface:

```bash
export CODONAMR_EMAIL="you@example.org"
codonamr run --ids ids.txt --out results/
```

where `ids.txt` is one NCBI accession per line:

```
NG_048025.1
NG_047831.1
NG_050430.1
```

That single command fetches the coding sequences, runs quality control and
reports what it rejected and why, identifies the source organism from the
records themselves, downloads the RefSeq reference or representative genome for
that organism, builds the host codon table and the CAI reference set from it,
computes every index in [section 5](#5-what-it-computes), draws the
length-matched null and the composition control, runs the permutation tests,
folds the mRNA if ViennaRNA is installed, validates the controls on this
dataset's own base composition, and writes `results/report.md` with the figures
and every intermediate table beside it.

Protein accessions work too. They are resolved to their coding nucleotide region
through the `/coded_by` qualifier of the protein record, which is exact.

From Python:

```python
from codonamr import analyse

result = analyse(ids=["NG_048025.1", "NG_047831.1"],
                 out_dir="results", email="you@example.org")
print(result["metrics"][["id", "enc", "gc3s", "cai_reference"]])
```

With no network at all, using your own sequences and your own host:

```bash
codonamr run --fasta my_genes.fasta --host host_cds.fasta --out results/
```

Indices only, no host comparison and no controls:

```bash
codonamr metrics --fasta my_genes.fasta --out metrics.tsv
```

---

## 5. What it computes

Every metric names the paper it implements, in the code and here. Where a
convention is ambiguous in the original, the choice made is stated in the
function's docstring, because different tools resolve those ambiguities
differently and the resulting numbers are not interchangeable. Formulae are in
[`docs/metrics.md`](docs/metrics.md).

| Metric | Function | Primary reference |
|---|---|---|
| Relative synonymous codon usage (RSCU) | `rscu` | Sharp PM, Tuohy TMF, Mosurski KR (1986) *Nucleic Acids Res* **14**:5125-5143 |
| Effective number of codons (ENC, Nc) | `enc` | Wright F (1990) *Gene* **87**:23-29 |
| Codon adaptation index (CAI) | `cai`, `build_w` | Sharp PM, Li WH (1987) *Nucleic Acids Res* **15**:1281-1295 |
| Frequency of optimal codons (Fop) | `fop`, `optimal_codons` | Ikemura T (1981) *J Mol Biol* **146**:1-21 |
| Synonymous codon usage order (SCUO) | `scuo` | Wan XF, Xu D, Kleinhofs A, Zhou J (2004) *BMC Evol Biol* **4**:19 |
| tRNA adaptation index (tAI) | `tai` | dos Reis M, Savva R, Wernisch L (2004) *Nucleic Acids Res* **32**:5036-5044 |
| Codon usage frequency similarity (CUFS) | `cufs` | Diament A, Pinter RY, Tuller T (2014) *Nat Commun* **5**:5876; metric from Endres DM, Schindelin JE (2003) *IEEE Trans Inf Theory* **49**:1858-1860 |
| GC, GC3s, GC12 | `gc_content`, `gc3s` | Wright F (1990) *Gene* **87**:23-29 for the GC3s convention |
| Neutrality plot (GC12 against GC3) | `report.plot_neutrality` | Sueoka N (1988) *Proc Natl Acad Sci USA* **85**:2653-2657 |
| PR2 bias plot | `report.plot_pr2` | Sueoka N (1995) *J Mol Evol* **40**:318-325 |
| Local mRNA folding energy | `dg_windowed` | Lorenz R et al. (2011) *Algorithms Mol Biol* **6**:26 (ViennaRNA); Hofacker IL et al. (1994) *Monatsh Chem* **125**:167-188 |
| Synonymous shuffle, codon order only | `randomise_shuffle` | Workman C, Krogh A (1999) *Nucleic Acids Res* **27**:4816-4822 |
| Dinucleotide-preserving synonymous shuffle | `randomise_shuffle_dinuc` | Katz L, Burge CB (2003) *Genome Res* **13**:2042-2051 |
| Empirical p-value correction | `empirical_z` | North BV, Curtis D, Sham PC (2002) *Am J Hum Genet* **71**:439-441 |
| tRNA gene copy numbers | `trna_copy_numbers_from_tRNAscan` | Chan PP, Lowe TM (2019) *Methods Mol Biol* **1962**:1-14 |
| Sequence retrieval and parsing | `fetch_cds`, `read_fasta` | Sayers EW et al. (2022) *Nucleic Acids Res* **50**:D20-D26; Cock PJA et al. (2009) *Bioinformatics* **25**:1422-1423 |

tAI needs tRNA gene copy numbers, which cannot be inferred from a coding
sequence. Supply tRNAscan-SE output for the host with `--trnascan`, or tAI is
reported as NA and the report says why.

---

## 6. What the controls do and why

This is the part of the package that is not a calculator.

### The length-matched null

ENC is estimated from the homozygosity of each synonymous family. When a family
holds only a few codons that estimate is biased upward and highly variable.
Below roughly 200 codons the bias is large enough to dominate any biological
signal, and below about 100 codons ENC is close to uninformative. The practical
consequence is that a short gene looks less biased than a long gene from the
same genome, with no difference in selection at all.

`length_matched_null` therefore compares each query gene against host genes of
its own length rather than against the genome mean. The bias then applies
equally to both sides and cancels. The pipeline reports a z score and an
empirical p-value against that matched distribution for ENC, SCUO, GC3s and CAI.
When too few host genes match on length, the tolerance is widened, and if genes
still have to be drawn with replacement, the output records it and the report
says those p-values are lower bounds. SCUO is bounded and behaves sensibly on
short sequences, so it is the safer index to lead with for short genes.

### The per nucleotide composition control

Controlling for GC alone treats A and T, and C and G, as interchangeable. They
are not. Replication-associated strand asymmetry, amino acid composition and
transcription-coupled repair all act on individual bases, and a gene acquired
from a compositionally different genome typically differs in more than its GC
total. `composition_control` runs a partial Spearman correlation on the A, C and
G fractions separately (T is dropped because the four sum to 1 and the design
matrix would otherwise be singular, and the remaining three span the same
space). It reports the uncontrolled correlation next to the controlled one. If
the partial correlation survives that and the GC-only one does not, it is the
GC-only result that was wrong.

The control removes monotone dependence on composition, not arbitrary
dependence. A U-shaped relationship with GC survives it. That limitation is
stated in the docstring rather than buried.

### The randomisation designs

Two nulls, holding different things fixed, because what a null preserves decides
what a significant result means.

`randomise_replace` changes codon composition while holding the peptide fixed,
so it answers what codon composition does. `randomise_shuffle` permutes codon
positions within each synonymous family, holding the peptide, the codon counts,
ENC, CAI and GC3s all exactly fixed, so anything it moves is attributable to
codon order alone. But a within-family permutation does not preserve
dinucleotide composition across codon junctions, and folding energy depends on
base stacking, so a positive result from that design can still be dinucleotide
composition rather than codon order. `randomise_shuffle_dinuc` adds the
dinucleotide constraint and **reports its residual**: a residual of zero means
the match is exact, and a non-zero residual means part of the effect could still
be composition. The report prints the residual rather than assuming it reached
zero.

---

## 7. Validating the implementations

The indices here are implemented from the primary literature rather than taken
from an established package, and the controls are the point of the package, so
both are tested rather than asserted.

### The controls are validated by injection

```bash
codonamr validate
```

This runs two experiments and prints the table.

`injection_recovery` plants a partial correlation of known size between two
variables that also share a composition-driven confounder, then analyses the
same data three ways: with no control, with a GC-only control, and with the per
nucleotide composition control. A control that works recovers the planted
effect. `false_positive_rate` plants nothing and reports how often an
association is found anyway. A control that works rejects at close to alpha.

A typical result at n = 200, a confounder of strength 0.8, and half of that
confounder carried by AT skew:

| experiment | control | mean r | true r | rejection rate |
|---|---|---|---|---|
| injection_recovery | no control | 0.56 | 0.30 | 1.00 |
| injection_recovery | GC only | 0.45 | 0.30 | 1.00 |
| injection_recovery | per-nucleotide composition | 0.30 | 0.30 | 0.99 |
| false_positive_rate | no control | 0.38 | 0.00 | 1.00 |
| false_positive_rate | GC only | 0.24 | 0.00 | 0.92 |
| false_positive_rate | per-nucleotide composition | 0.03 | 0.00 | 0.07 |

Read the second block first. With no control, a completely null association is
declared significant every single time. With a GC-only control it is still
declared significant 92 times in 100. With the per nucleotide control the false
positive rate falls to roughly the nominal 5 per cent, and the first block shows
that this is not bought by destroying real signal: the planted effect of 0.30
comes back as 0.30.

Set `--skew 0` and the GC-only control performs as well as the per nucleotide
one, which is the honest statement of when the extra covariates matter: they
matter when the confounding is strand-asymmetric, and not otherwise. The
pipeline runs both experiments on your own dataset's base composition and writes
the result to `validation.tsv` on every run.

### The indices are tested against their exact invariants

```bash
python -m pytest tests/ -q
```

The suite asserts properties that are true by definition of the published
method, not that the output has not changed since yesterday:

- ENC is exactly 20 for a sequence using one codon per amino acid, and in the
  top of its range for uniform synonymous usage;
- RSCU sums to the family size within every synonymous family, and returns
  `None`, not zero, for a family that was never used;
- `randomise_shuffle` preserves codon counts, ENC, CAI and GC3s exactly, and
  preserves the encoded peptide;
- `randomise_shuffle_dinuc` preserves codon counts and its reported residual
  equals the true L1 distance from the native junction profile;
- CAI of a maximally biased reference set against itself is exactly 1;
- `spearman_partial` recovers a planted partial correlation to within the known
  attenuation of the rank transform;
- every QC failure mode returns its documented reason string;
- the tAI wobble weights equal the published constants of dos Reis et al. (2004)
  exactly;
- the same seed gives byte-identical output tables, on one thread or several.

Tests that need the network are marked `network` and are deselected by default,
so the suite runs offline. Run them with `pytest -m network` and
`CODONAMR_EMAIL` set.

If you intend to rely on these numbers, cross-check a subset against CodonW,
EMBOSS `cusp` and `chips`, or CAIcal. Discrepancies are worth reporting as
issues, and [section 12](#12-support) says where.

---

## 8. Command line

| command | what it does |
|---|---|
| `codonamr run` | the full analysis: fetch, QC, metrics, host comparison, controls, randomisation, report |
| `codonamr metrics` | per-sequence indices only, no host comparison and no network |
| `codonamr validate` | the control validation of [section 7](#7-validating-the-implementations) |
| `codonamr version` | print the version |

Useful options for `run`:

| option | effect |
|---|---|
| `--fasta FILE` | local sequences instead of, or alongside, `--ids` |
| `--host TAXID_OR_NAME_OR_FASTA` | override the host inferred from the records |
| `--trnascan FILE` | tRNAscan-SE output for the host, which enables tAI |
| `--no-folding` | skip mRNA folding even if ViennaRNA is installed |
| `--permutations N` | permutations per gene, default 100 |
| `--threads N` | worker processes, `1` runs everything in this process |
| `--seed N` | master seed, default 0 |
| `--temperature C` | folding temperature in Celsius |

`codonamr run --help` lists the rest.

---

## 9. Output files

Every intermediate is a tab-separated file, so no result is trapped inside the
program.

| file | contents |
|---|---|
| `report.md` | the readable summary, with the QC counts, the missing-value table, the controls and the figures |
| `qc.tsv`, `qc_summary.tsv` | every sequence with its QC verdict, and the counts per rejection reason |
| `metrics.tsv` | per-gene indices |
| `per_gene_all.tsv` | metrics joined to the randomisation and folding results |
| `rscu_per_gene.tsv`, `codon_counts.tsv`, `aa_composition.tsv` | the raw tables behind the indices |
| `host_codon_table.tsv` | host counts, RSCU, genome-wide w, reference-set w and the optimal codon flags |
| `host_gene_metrics.tsv` | the host's own genes, the background distribution the query is read against |
| `controls_length_matched.tsv` | per gene and per metric: observed, null mean, z, empirical p, pool size, and whether resampling was needed |
| `controls_composition.tsv` | raw and partial correlations for each tested pair |
| `randomisation.tsv` | permutation statistics, the dinucleotide residual and the folding null |
| `validation.tsv` | the injection recovery and false positive rate tables for this dataset |
| `manifest.json` | the full run record, including the seed, the host assembly and every note |
| `figures/` | 300 dpi PNG and vector PDF, in a colourblind-safe palette |

Figures: ENC against GC3s with the Wright expected curve, the neutrality plot,
the PR2 bias plot, an RSCU heatmap of the genes against the host, and a
gene-against-host compatibility panel.

The report states the QC rejection counts and lists every metric that returned
`None`, with the reason. Nothing is dropped silently.

---

## 10. Limitations

These are properties of the methods, not bugs, and none of them are fixable by
better code.

**Deposition bias is not prevalence.** The number of accessions for a resistance
gene reflects what was sequenced, funded, and deposited. Counting records is not
an epidemiological estimate, and no output of this package should be read as
one.

**CAI reference sets are not comparable across genomes.** CAI is defined against
a set of highly expressed genes from one genome. A CAI computed against
*Escherichia coli* ribosomal proteins and one computed against *Acinetobacter
baumannii* ribosomal proteins are not on the same scale, and the difference
between them is not a biological quantity. Scoring a gene against its own
genome's reference set is also partly circular, since the gene contributed to
the genome that defined the set. Both caveats bite hardest for exactly the case
this package is built for, a mobile gene compared across hosts.

**ENC on short genes.** ENC is biased upward and noisy below roughly 200 codons
and close to uninformative below about 100. The length-matched null exists
because of this, but a control makes the comparison fair, it does not create
information that the sequence does not contain. For short genes, prefer SCUO and
read ENC with the reported z score rather than on its own.

**Folding at non-physiological temperatures.** Folding energies come from a
nearest-neighbour thermodynamic model applied to naked RNA, with no ribosomes,
no proteins and no ions beyond the model's assumptions. They are a comparative
statistic between sequences, not a physical energy. Be especially careful with
thermophiles: at 80 degrees Celsius the model retains almost no structure, so
the statistic loses most of its resolution, and comparing a thermophile folded
at 80 with a mesophile folded at 37 compares two different things.

**Host resolution can be approximate.** If no reference or representative genome
exists for the exact taxid, the match falls back to the species and then to the
genus, and the codon table then comes from a relative rather than the strain.
This is recorded in `manifest.json` and stated in the report, but a comparison
against the wrong host is not evidence of anything.

**Gene copy number is a proxy for tRNA abundance.** It is a good proxy in fast
growing bacteria and a poor one in slow growers and in eukaryotes. The tAI
wobble constants were fitted to two organisms, *E. coli* and *S. cerevisiae*,
and are used unchanged everywhere, which is standard practice and still an
assumption.

---

## 11. Documentation

- [`docs/index.md`](docs/index.md), the landing page
- [`docs/metrics.md`](docs/metrics.md), the formula and primary reference for
  every metric
- [`docs/api.md`](docs/api.md), the public API
- [`CONTRIBUTING.md`](CONTRIBUTING.md), how to report a bug or propose a change

---

## 12. Support

Questions, bug reports and feature requests go to the GitHub issue tracker:

<https://github.com/bioforkgithub/codonamr/issues>

A useful bug report includes the command or the Python call you ran, the
`manifest.json` from the run if there is one, the versions of codonamr and
Python, and what you expected instead. Numerical discrepancies against CodonW,
EMBOSS or CAIcal are welcome and will be treated as bugs until shown otherwise.

Please do not paste unpublished sequence data into an issue. A minimal synthetic
sequence that reproduces the problem is better for both of us.

---

## 13. Citation

If you use this software, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff), and cite the methods it implements. The primary
references are listed in [section 5](#5-what-it-computes) and given with their
formulae in [`docs/metrics.md`](docs/metrics.md).

---

## 14. License

Released under the **MIT License**, see [`LICENSE`](LICENSE).
Copyright (c) 2026 Manish Prakash Victor.

You may use, modify and redistribute this code freely, including commercially,
provided the copyright notice and licence text are retained. The code is
provided without warranty.

Note that the licence covers this software only. The sequence records it
downloads are subject to the terms of NCBI and RefSeq.
