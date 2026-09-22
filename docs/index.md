# codonamr documentation

Codon usage analysis of antimicrobial resistance genes, with a length-matched
null and a per nucleotide composition control that are themselves validated.

You supply sequence identifiers. The package fetches the coding sequences,
resolves the source organism to a RefSeq reference genome, builds the host codon
table and CAI reference set from it, computes the standard indices, applies the
controls, runs the permutation tests, and writes a Markdown report with figures
and every intermediate as a TSV.

## Contents

- [Metrics](metrics.md): the formula and primary reference for every metric.
- [API](api.md): the public functions, their arguments and what they return.
- [README](../README.md): statement of need, how this differs from existing
  tools, what the controls do and why, and the standing limitations.
- [Contributing](../CONTRIBUTING.md): reporting a bug, proposing a change,
  running the tests.

## Install

```bash
git clone https://github.com/bioforkgithub/codonamr
cd codonamr
pip install -e ".[dev]"
```

`pip install codonamr` will work once the first release is on PyPI.

Python 3.9 or newer. Requires biopython, numpy, scipy, pandas and matplotlib.
Local mRNA folding is optional because ViennaRNA is a large install:

```bash
pip install "codonamr[folding]"
```

NCBI asks that every request carries a contact address. Set it once:

```bash
export CODONAMR_EMAIL="you@example.org"
```

Downloads are cached on disk under the accession, so a repeat run makes no
network calls. Set `CODONAMR_CACHE` to choose where, and `CODONAMR_API_KEY` to
use an NCBI API key, which raises the rate limit from 3 to 10 requests per
second.

## Quickstart

Everything from identifiers:

```bash
codonamr run --ids ids.txt --out results/
```

`ids.txt` holds one NCBI accession per line. Protein accessions are accepted and
are resolved to their coding nucleotide region through the `/coded_by`
qualifier. The run writes `results/report.md`, the figures, and a TSV for every
intermediate.

From Python:

```python
from codonamr import analyse

result = analyse(ids=["NG_048025.1"], out_dir="results",
                 email="you@example.org")
result["metrics"][["id", "n_codons", "enc", "gc3s", "cai_reference"]]
```

Offline, with your own sequences and your own host:

```bash
codonamr run --fasta my_genes.fasta --host host_cds.fasta --out results/
```

Indices only:

```bash
codonamr metrics --fasta my_genes.fasta --out metrics.tsv
```

Check that the controls do what they claim:

```bash
codonamr validate
```

## What a run produces

`report.md` states the QC rejection counts and every metric that returned
`None`, with the reason, so nothing is dropped silently. Beside it:
`metrics.tsv`, `per_gene_all.tsv`, `qc.tsv`, `qc_summary.tsv`,
`host_codon_table.tsv`, `controls_length_matched.tsv`,
`controls_composition.tsv`, `randomisation.tsv`, `validation.tsv`,
`manifest.json`, and `figures/` with 300 dpi PNG and vector PDF versions of the
ENC against GC3s plot, the neutrality plot, the PR2 plot, the RSCU heatmap and
the gene-against-host compatibility panel.

Full descriptions are in section 9 of the [README](../README.md).

## Reproducibility

The same seed gives the same output. Per-gene random seeds are derived from the
master seed and the record identifier rather than from call order, so the result
does not depend on how many worker processes were used or on the order of the
input. `--threads 1` runs everything in one process. `manifest.json` records the
seed, the host assembly that was chosen, how it was matched, and every note the
run generated.
