# Contributing to codonamr

Contributions are welcome, including bug reports, corrections to the
implementations, new metrics with their primary references, and documentation
fixes. This file says how each of those works.

The package implements published methods in house rather than wrapping
established binaries, so **a report that a number disagrees with CodonW, EMBOSS
or CAIcal is one of the most useful things you can send**. Those are treated as
bugs until shown otherwise.

---

## 1. Reporting a bug

Open an issue at
<https://github.com/bioforkgithub/codonamr/issues>.

A useful report has:

- what you ran, the exact command line or the Python call;
- what you expected and what happened instead;
- the `manifest.json` from the run, if there was one, since it records the seed,
  the host assembly that was chosen and every note the run generated;
- versions: `codonamr version`, `python --version`, and the versions of numpy,
  scipy, pandas and biopython;
- a minimal example that reproduces the problem.

Please do not paste unpublished sequence data into an issue. A short synthetic
sequence that triggers the same behaviour is better, and `tests/conftest.py` has
generators for making one.

For a numerical discrepancy, say which tool you compared against, which version
of it, and give both numbers. Different tools resolve the ambiguities in the
original papers differently, so the first question is always whether the two
programs are computing the same quantity. The conventions this package follows
are stated in each function's docstring and in `docs/metrics.md`.

---

## 2. Proposing a change

Open an issue before starting work on anything larger than a typo. That is not
a formality: it avoids two people writing the same patch, and for a new metric
it is where the choice of primary reference gets settled.

Then:

1. Fork the repository and branch from `main`.
2. Make the change, with tests.
3. Run the test suite, see section 4.
4. Open a pull request describing what changed and why. Link the issue.

A pull request that changes a computed value must say which value changed, by
how much, and why the new one is right. Include the reference if the answer is
in a paper.

### Adding a metric

A new metric needs all of these:

- the primary reference, named in the module docstring and in the function
  docstring, and added to the table in `docs/metrics.md` with the formula as
  coded;
- a statement in the docstring of any convention the original leaves ambiguous,
  and which way this implementation resolves it;
- `None` rather than a placeholder number when the metric cannot be computed,
  plus an entry in `_MISSING_REASONS` in `report.py` so the report can say why;
- at least one test asserting an exact invariant of the method, not a smoke
  test. "ENC is exactly 20 when one codon is used per amino acid" is an
  invariant. "It runs without raising" is not.

---

## 3. Coding conventions

- **Python 3.9 and newer.** No syntax that needs a later version.
- **Dependencies.** The standard library plus biopython, numpy, scipy, pandas
  and matplotlib. ViennaRNA stays an optional extra, imported lazily, and the
  codon usage half of the package must keep working without it. Adding a new
  required dependency needs a discussion in an issue first.
- **No em-dashes** anywhere in code, comments, docstrings or documentation. Use
  a comma, a parenthesis or a full stop.
- **No hardcoded paths, no credentials, no private data.** This repository is
  public and contains code only. Cache and output locations come from arguments
  or environment variables.
- **Docstrings** are numpydoc-style. State caveats in prose where a choice is
  ambiguous. The docstring is where a reader finds out that CAI values from two
  genomes are not comparable, so write it there rather than assuming it is
  common knowledge.
- **Every network call is polite.** An email address is mandatory, requests are
  throttled to the NCBI limit, failures are retried with backoff, and responses
  are cached on disk so a rerun does not refetch.
- **Determinism.** Anything stochastic takes a seed or a `random.Random`.
  Per-item seeds are derived from the master seed and the item's identifier, not
  from call order, so results do not depend on the number of worker processes.
  `hash()` is salted per process in Python 3 and must not be used for this.
- **Write every intermediate as TSV.** No result should be reachable only from
  inside the program.
- Four-space indentation, lines under 80 characters where that does not hurt
  readability, and imports grouped standard library, third party, local.

There is no enforced formatter. Match the surrounding style.

---

## 4. Running the tests

```bash
pip install -e ".[dev]"
python -m pytest tests/ -q
```

The suite runs offline. Tests that need the network are marked `network` and are
deselected by default through `addopts` in `pyproject.toml`. To run them:

```bash
export CODONAMR_EMAIL="you@example.org"
python -m pytest tests/ -m network
```

Coverage, if you want it:

```bash
python -m pytest tests/ --cov=codonamr --cov-report=term-missing
```

The same suite runs in GitHub Actions on Python 3.9, 3.11 and 3.12 for every
push and pull request, offline. A pull request is expected to be green there.

Tests should assert invariants of the method. The existing suite is the model:
ENC is exactly 20 for one codon per amino acid, RSCU sums to the family size,
`randomise_shuffle` preserves codon counts, ENC, CAI and GC3s exactly, and the
tAI wobble weights equal the published constants.

---

## 5. How a contribution gets reviewed

Review is by the maintainer. What is checked, in order:

1. **Is it correct?** For anything touching a metric or a control, this means
   checking the implementation against the cited paper, not just reading the
   diff.
2. **Is there a test that would fail without the change?** A bug fix needs a
   regression test. A new metric needs an invariant test.
3. **Does the test suite pass** on the supported Python versions, offline.
4. **Does the documentation match?** A changed formula means a changed
   `docs/metrics.md`. A changed signature means a changed `docs/api.md`.
5. **Are the caveats stated?** If the change introduces a quantity that is not
   comparable across genomes, or that is unstable on short sequences, the
   docstring has to say so.
6. **Conventions**, from section 3, including the em-dash rule.

Expect questions about the choice of reference and about what a new statistic
holds fixed. Those are not obstruction: what a null model preserves decides what
a significant result means, and getting that wrong is the failure mode this
package exists to avoid.

Small, focused pull requests are reviewed faster than large ones. If a change
turns out to be larger than expected, split it.

---

## 6. Code of conduct

Participation in this project is governed by
[`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md), which is the Contributor Covenant
version 2.1.

---

## 7. Licence

By contributing you agree that your contribution is licensed under the MIT
Licence, the same terms as the rest of the project. See [`LICENSE`](LICENSE).
