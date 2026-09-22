# Metrics: formulae and primary references

Every metric implemented in `codonamr`, with the formula as it is coded and the
paper it comes from. Where the original leaves a convention undecided, the
choice made here is stated, because different tools resolve those ambiguities
differently and the resulting numbers are not interchangeable.

Notation used throughout: a synonymous family is the set of codons encoding one
amino acid, `k` is the family size, `n_i` is the count of codon `i`, and `N` is
the number of codons observed in the family. Stop codons are excluded
everywhere. Methionine and tryptophan have one codon each and carry no
synonymous information, so they are excluded from every index that measures
synonymous choice.

---

## Relative synonymous codon usage (RSCU)

Sharp PM, Tuohy TMF, Mosurski KR (1986) *Nucleic Acids Res* **14**:5125-5143.

```
RSCU_i = n_i / ( (1/k) * sum_j n_j )
       = k * n_i / N
```

The observed count divided by the count expected if every codon in the family
were used equally. RSCU = 1 means no bias. Within a family the values sum to
`k` by construction, which is the invariant the test suite checks.

Single-codon families are 1 by definition. A family with no observed codons
returns `None` for every member rather than 0, so that "not seen" stays
distinguishable from "seen and not used". Function: `codonamr.rscu`.

---

## Effective number of codons (ENC, Nc)

Wright F (1990) *Gene* **87**:23-29.

For each family, the homozygosity estimate is

```
F_hat = ( N * sum_i (n_i / N)^2  -  1 ) / ( N - 1 )
```

and

```
Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
```

where `F_k` is the mean of `F_hat` over the families with degeneracy `k`. The
constant 2 is the contribution of Met and Trp, one codon each.

Conventions followed here, matching CodonW:

- families with `N < 2` carry no information and are dropped;
- if the single threefold family (isoleucine) is unusable, `F3` is interpolated
  as the mean of `F2` and `F4`;
- a degeneracy class with no usable family contributes nothing, and its term is
  dropped from the sum rather than treated as zero;
- the result is capped at 61.

Nc ranges from 20, one codon used per amino acid, to 61, all synonymous codons
used equally. The estimator slightly exceeds 61 at finite sequence length, which
is why it is capped.

`enc` returns `None` below 20 codon-equivalents of information rather than a
precise-looking number computed from almost nothing. Function: `codonamr.enc`.

**Expected value under no selection** (Wright 1990), plotted as the curve on the
ENC against GC3s figure, with `s` the GC3s value:

```
Nc_expected = 2 + s + 29 / ( s^2 + (1 - s)^2 )
```

Function: `codonamr.report.wright_expected_enc`.

---

## Codon adaptation index (CAI)

Sharp PM, Li WH (1987) *Nucleic Acids Res* **15**:1281-1295.

Relative adaptiveness, from a reference set of highly expressed genes:

```
w_i = f_i / f_max        within each synonymous family
```

where `f_max` is the count of the commonest codon in that family in the
reference set. Codons never seen in the reference set are given `0.5 / f_max`
rather than zero, the standard correction that keeps the geometric mean finite.

```
CAI = ( prod_over_codons w_i ) ^ (1 / L)
```

the geometric mean of `w` over the `L` informative codons of the gene, which
excludes Met, Trp and stops.

The reference set defines what CAI means. Values computed against different
genomes are not on the same scale. Functions: `codonamr.build_w`,
`codonamr.cai`, `codonamr.reference_set`.

---

## Frequency of optimal codons (Fop)

Ikemura T (1981) *J Mol Biol* **146**:1-21.

```
Fop = (number of optimal codons) / (number of codons in degenerate families)
```

Ikemura defined optimal codons against measured tRNA abundances. Without
expression data the operational definition used here is RSCU > 1 in the host
genome, which is a proxy for that and comes apart from it in genomes with strong
mutational GC skew. Functions: `codonamr.fop`, `codonamr.optimal_codons`.

---

## Synonymous codon usage order (SCUO)

Wan XF, Xu D, Kleinhofs A, Zhou J (2004) *BMC Evol Biol* **4**:19.

Per family, with `p_i = n_i / N`:

```
H      = - sum_i p_i * ln(p_i)
H_max  = ln(k)
O      = (H_max - H) / H_max
```

and the gene value is the mean of `O` over families, weighted by the number of
codons in each family:

```
SCUO = sum_families (O * N) / sum_families N
```

0 means codons within each family are used evenly, 1 means a single codon is
used throughout. Unlike ENC it is bounded and behaves sensibly on short
sequences, which makes it the safer index for short resistance genes. Function:
`codonamr.scuo`.

---

## tRNA adaptation index (tAI)

dos Reis M, Savva R, Wernisch L (2004) *Nucleic Acids Res* **32**:5036-5044.

Absolute adaptiveness of codon `i`:

```
W_i = sum_j ( 1 - s_ij ) * tGCN_j
```

summed over the tRNAs `j` that can decode codon `i`, where `tGCN_j` is the gene
copy number of tRNA `j` and `s_ij` is the selective constraint on that pairing.
The constants, with the anticodon wobble base 34 first and the codon base 3
second:

| pairing | s | meaning |
|---|---|---|
| I:U | 0.00 | A34, read as inosine, on codon U3 |
| G:C | 0.00 | Watson-Crick |
| U:A | 0.00 | Watson-Crick |
| C:G | 0.00 | Watson-Crick |
| G:U | 0.41 | wobble |
| I:C | 0.28 | wobble |
| I:A | 0.9999 | possible in principle, effectively never used |
| U:G | 0.68 | wobble |
| L:A | 0.89 | lysidine-modified C34 of bacterial tRNA-Ile2 on codon AUA |

Weights are then normalised and the zeros are filled, which is what keeps the
geometric mean finite:

```
w_i = W_i / max(W)
w_i = geometric mean of the non-zero w   if W_i = 0

tAI = ( prod_over_codons w_i ) ^ (1 / L)
```

Met and Trp are excluded by default, as in the original, because they carry no
synonymous choice. Including them shifts the index by a genome-dependent
constant, so values computed with and without are not comparable.

The bacterial lysidine tRNA-Ile shares the anticodon string CAU with tRNA-Met,
so `trna_copy_numbers_from_tRNAscan` stores it under the separate key
`CAT_Ile2`. Functions: `codonamr.tai`,
`codonamr.trna_copy_numbers_from_tRNAscan`.

The `s` values were fitted to *E. coli* and *S. cerevisiae*. They are empirical
constants from two organisms, not physical ones.

---

## Codon usage frequency similarity (CUFS)

Diament A, Pinter RY, Tuller T (2014) *Nat Commun* **5**:5876, using the metric
of Endres DM, Schindelin JE (2003) *IEEE Trans Inf Theory* **49**:1858-1860.

RSCU is a per-family quantity, not a probability distribution, so each family is
first converted back to within-family frequencies and the families are then
weighted equally over the `m` families both tables define:

```
p_i = RSCU_i / (k * m)
```

Then, with `M = (P + Q) / 2` and base-2 logarithms,

```
JSD(P, Q) = 0.5 * sum_i p_i log2(p_i / m_i) + 0.5 * sum_i q_i log2(q_i / m_i)
CUFS      = sqrt( JSD(P, Q) )
```

0 means identical usage, 1 means the two never choose the same codon. Despite
the name it is a distance, not a similarity. Equal weighting of families is a
choice: weighting by amino-acid usage instead would make the result partly a
measure of protein composition. Function: `codonamr.cufs`.

`codonamr.rscu_distance` offers the plain Euclidean, Manhattan, cosine and
correlation distances between two RSCU tables for the same comparison.

---

## Composition measures

**GC content.** The fraction of G and C over the whole sequence.

**GC3s.** GC at synonymous third positions. Third positions of Met, Trp and stop
codons are excluded, so this is GC3s and not plain GC3. Published GC3 values
should be checked for which was meant.

**GC12 and GC3, for the neutrality plot.** Sueoka N (1988) *Proc Natl Acad Sci
USA* **85**:2653-2657. Here all sense codons are used, including the
single-codon families, because the neutrality plot is a statement about
mutational pressure across the gene rather than about synonymous choice. The
regression slope of GC12 on GC3 estimates the share of variation attributable to
genome-wide mutational pressure: a slope near 1 means mutation dominates, a
slope near 0 means selection is acting on the coding positions.

**PR2 bias.** Sueoka N (1995) *J Mol Evol* **40**:318-325. At fourfold
degenerate third positions only:

```
x = A3 / (A3 + T3)
y = G3 / (G3 + C3)
```

Under no bias both coordinates are 0.5. Distance from that centre measures
strand asymmetry. A gene can sit at GC3 = 0.5 and still be far off centre here,
which is the concrete reason a GC-only control is not sufficient.

Functions: `codonamr.gc_content`, `codonamr.gc3s`, and the plotting functions in
`codonamr.report`.

---

## Local mRNA folding energy

Lorenz R et al. (2011) *Algorithms Mol Biol* **6**:26 (ViennaRNA); Hofacker IL
et al. (1994) *Monatsh Chem* **125**:167-188.

```
dG_mean = mean over windows of MFE( sequence[i : i + window] )
```

sliding by `step`. The defaults are a 40-nucleotide window with a step of 3.
ViennaRNA is an optional dependency, and `codonamr.folding.available()` reports
whether it is installed. Function: `codonamr.dg_windowed`.

Folding energies are a comparative statistic between sequences, not a physical
energy: the model assumes naked RNA with no ribosomes and no proteins.

---

## Randomisation designs

**Design A, synonymous replacement.** Each codon is replaced by a synonymous
one, with the family's preferred codon drawn with marginal probability `p` and
the rest spread evenly. `p = 1` gives maximal bias, ENC near 20, and `p = 1/k`
gives uniform usage, ENC near 61. Holds the peptide fixed and changes codon
composition, so it answers what codon composition does. Function:
`codonamr.randomise_replace`.

**Design B, within-family permutation.** Workman C, Krogh A (1999) *Nucleic
Acids Res* **27**:4816-4822. Codon positions are permuted within each synonymous
family. The peptide, the codon counts, ENC, CAI and GC3s are all invariant by
construction, so anything that moves is attributable to codon order. Function:
`codonamr.randomise_shuffle`.

**Design B with dinucleotide control.** Katz L, Burge CB (2003) *Genome Res*
**13**:2042-2051. Design B does not preserve dinucleotide composition across
codon junctions, and folding energy depends on base stacking, so a result from
design B alone can reflect dinucleotide composition rather than codon order.
This design starts from an unconstrained permutation and accepts only synonymous
swaps that reduce the L1 distance between the shuffled and native
junction-dinucleotide profiles. It returns the sequence and the residual L1
distance. A residual of 0 means the composition is matched exactly. The residual
is reported, never assumed. Function: `codonamr.randomise_shuffle_dinuc`.

Codon-internal dinucleotides are fixed by codon composition alone, so they are
invariant under any within-family permutation. Junction dinucleotides, base 3 of
codon `i` with base 1 of codon `i+1`, are the only part of dinucleotide
composition that codon order can change. Function: `codonamr.junction_counts`.

---

## Controls

**Length-matched null.** Host genes are sampled to match the query gene's codon
length, and the query's index is standardised against that matched
distribution rather than against the genome mean. The tolerance starts at 20 per
cent and widens if too few genes match. ENC's small-sample bias is documented in
Fuglsang A (2004) *Biochem Biophys Res Commun* **317**:957-964. Function:
`codonamr.length_matched_null`.

**Empirical p-value.** North BV, Curtis D, Sham PC (2002) *Am J Hum Genet*
**71**:439-441.

```
p = (r + 1) / (m + 1)
```

with `r` the number of null values at least as extreme as the observation and
`m` the number of draws. It can never be exactly zero, so with 100 permutations
the smallest attainable two-sided p is about 0.02. That is a property of the
design, not of the data. Function: `codonamr.controls.empirical_z`.

**Partial Spearman correlation.** Kendall MG, Stuart A (1973) *The Advanced
Theory of Statistics*, volume 2. Everything is rank-transformed, `x` and `y` are
regressed on the ranked covariates, and Pearson's correlation is taken between
the residuals:

```
r_partial = corr( resid(rank(x) ~ rank(C)), resid(rank(y) ~ rank(C)) )
t         = r * sqrt( df / (1 - r^2) ),  df = n - rank(design) - 1
```

The covariates are the A, C and G fractions of each sequence. T is dropped
because the four fractions sum to 1 and the design matrix would otherwise be
singular; the remaining three span the same space, so nothing is lost. This
removes monotone dependence on composition, not arbitrary dependence: a U-shaped
relationship with GC survives it. Functions: `codonamr.spearman_partial`,
`codonamr.composition_control`.

**Validation of the controls.** `codonamr.injection_recovery` plants a partial
correlation of known size in data with a composition confounder, part of which
is carried by AT skew and is therefore invisible to a GC-only control, and
reports what each control returns. `codonamr.false_positive_rate` plants nothing
and reports how often an association is found anyway. See section 7 of the
README for a worked table.

---

## Quality control

A sequence is rejected with one of these reasons, all reported and none applied
silently: `too_short`, `not_multiple_of_three`, `ambiguous_bases`,
`no_stop_codon`, `internal_stop`, `bad_start`. A passing sequence is returned as
uppercase DNA with the terminal stop removed. Function: `codonamr.check_cds`.
