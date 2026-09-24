"""Codon-usage indices, implemented from the primary literature.

Each function names the paper it implements. Where a convention is ambiguous in
the original, the choice made here is stated in the docstring, because different
tools resolve those ambiguities differently and the resulting numbers are not
interchangeable.

References
----------
ENC   Wright F (1990) Gene 87:23-29.
CAI   Sharp PM, Li WH (1987) Nucleic Acids Res 15:1281-1295.
RSCU  Sharp PM, Tuohy TMF, Mosurski KR (1986) Nucleic Acids Res 14:5125-5143.
SCUO  Wan XF, Xu D, Kleinhofs A, Zhou J (2004) BMC Evol Biol 4:19.
Fop   Ikemura T (1981) J Mol Biol 146:1-21.
"""
import math
from collections import defaultdict

from .genetic_code import CODON2AA, CODONS, DEGENERATE, FAMILY, codon_list

__all__ = ["codon_counts", "rscu", "enc", "build_w", "cai", "gc_content",
           "gc3s", "fop", "scuo", "aa_composition"]


def codon_counts(seq, drop_stops=True):
    """Counts of every sense codon, including zeros."""
    counts = {c: 0 for c in CODONS}
    for c in codon_list(seq):
        if c in counts:
            counts[c] += 1
    if drop_stops:
        counts = {c: n for c, n in counts.items() if CODON2AA[c] != "*"}
    return counts


def aa_composition(seq):
    """Amino-acid counts, stops excluded."""
    out = defaultdict(int)
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa and aa != "*":
            out[aa] += 1
    return dict(out)


# ------------------------------------------------------------------- RSCU
def rscu(seq_or_counts):
    """Relative synonymous codon usage (Sharp, Tuohy & Mosurski 1986).

    RSCU_ij = observed count / (count expected if all codons in the family were
    used equally). A value of 1 means no bias. Single-codon families are always
    1 by definition and are returned as such.

    Families with no observed codons return ``None`` for every member, rather
    than 0, so that "not seen" is distinguishable from "seen and not used".
    """
    counts = seq_or_counts if isinstance(seq_or_counts, dict) else codon_counts(seq_or_counts)
    out = {}
    for aa, fam in FAMILY.items():
        total = sum(counts.get(c, 0) for c in fam)
        k = len(fam)
        for c in fam:
            out[c] = None if total == 0 else counts.get(c, 0) * k / total
    return out


# -------------------------------------------------------------------- ENC
def enc(seq, min_expected=20):
    """Effective number of codons, Wright (1990).

    F-hat for a synonymous family of size n is ``(N*sum(p_i^2) - 1)/(N - 1)``,
    where N is the number of codons observed in that family. Families with
    N < 2 carry no information and are dropped.

    ``Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6``, with F_k the mean of F-hat over the
    families of degeneracy k.

    Conventions followed here, matching CodonW:

    * If the single 3-fold family (Ile) is unusable, F3 is interpolated as the
      mean of F2 and F4.
    * A degeneracy class with no usable family contributes nothing, and its
      term is dropped from the sum rather than treated as zero.
    * The result is capped at 61.

    ``min_expected`` is compared against the number of SYNONYMOUS FAMILIES that
    contributed, not against a codon count: the counter reaches 2 + 9 + 1 + 5 + 3
    = 20 when all four degeneracy classes are usable, so the default of 20 makes
    ``None`` mean "at least one whole degeneracy class was unusable". Sequence
    length never enters this function. A 36-codon input returns a number.

    This is therefore NOT a length or stability guard, and a non-``None`` return
    is not evidence that ENC is stable for that sequence. ENC destabilises when
    individual synonymous families carry few codons, which happens long before a
    whole class disappears; callers that need a length criterion (the usual rule
    of thumb is about 200 codons) must apply it themselves, and callers comparing
    ENC between sequences should note that the estimator carries a sampling sd of
    roughly 1.6 units and a composition-dependent downward offset of 1.5 to 4.3
    units even at 540 codons (measured by parametric bootstrap on the mcr set,
    ``classes/01_colistin_mcr/scripts/102_referee_enc_stability_and_composition.py``).

    Corrected 2026-09-24: this docstring previously said ``None`` was returned
    "when fewer than ``min_expected`` codon-equivalents of information are
    available ... Short genes hit this often". That described behaviour the code
    does not implement. The code is unchanged; only the description is.
    """
    counts = defaultdict(lambda: defaultdict(int))
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa and aa != "*":
            counts[aa][c] += 1

    byclass = defaultdict(list)
    for aa, fam in DEGENERATE.items():
        obs = counts.get(aa)
        if not obs:
            continue
        n = sum(obs.values())
        if n < 2:
            continue
        s = sum((v / n) ** 2 for v in obs.values())
        f = (n * s - 1.0) / (n - 1.0)
        if f <= 0:                       # even usage at small n can give f<=0
            continue
        byclass[len(fam)].append(f)

    m = {k: sum(v) / len(v) for k, v in byclass.items()}
    if 3 not in m and 2 in m and 4 in m:
        m[3] = (m[2] + m[4]) / 2.0

    total = expected = 2.0               # Met and Trp each contribute exactly 1
    for deg, n_families in ((2, 9), (3, 1), (4, 5), (6, 3)):
        if m.get(deg, 0) > 0:
            total += n_families / m[deg]
            expected += n_families
    if expected < min_expected:
        return None
    return min(total, 61.0)


# -------------------------------------------------------------------- CAI
def build_w(ref_seqs):
    """Relative adaptiveness ``w`` from a reference set of highly expressed genes.

    ``w_i = f_i / f_max`` within each synonymous family (Sharp & Li 1987).
    Codons never seen in the reference set are given ``0.5/f_max`` rather than
    zero, which is the standard correction to keep the geometric mean finite.

    The reference set defines what CAI means. A CAI computed against one genome's
    ribosomal proteins is not comparable to one computed against another's, and
    scoring a gene against its own genome's reference set is partly circular.
    Both caveats matter when comparing mobile genes across hosts.
    """
    counts = defaultdict(lambda: defaultdict(int))
    for s in ref_seqs:
        for c in codon_list(s):
            aa = CODON2AA.get(c)
            if aa and aa != "*":
                counts[aa][c] += 1
    w = {}
    for aa, fam in FAMILY.items():
        if len(fam) == 1:
            w[fam[0]] = 1.0
            continue
        obs = {c: counts[aa].get(c, 0) for c in fam}
        mx = max(obs.values())
        for c in fam:
            w[c] = (obs[c] / mx) if obs[c] > 0 else (0.5 / mx if mx else 0.01)
    return w


def cai(seq, w):
    """Codon adaptation index: geometric mean of ``w`` over informative codons.

    Met, Trp and stop codons are excluded, since they carry no synonymous choice.
    """
    total, n = 0.0, 0
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa is None or aa == "*" or len(FAMILY[aa]) == 1:
            continue
        wi = w.get(c, 0.01)
        total += math.log(wi if wi > 0 else 0.01)
        n += 1
    return math.exp(total / n) if n else None


# ------------------------------------------------------------ composition
def gc_content(seq):
    """Overall GC fraction."""
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else None


def gc3s(seq):
    """GC at synonymous third positions.

    Third positions of Met, Trp and stop codons are excluded, so this is GC3s
    and not plain GC3. Comparisons with published GC3 values should check which
    was meant.
    """
    g = n = 0
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa is None or aa == "*" or len(FAMILY[aa]) == 1:
            continue
        n += 1
        if c[2] in "GC":
            g += 1
    return g / n if n else None


def fop(seq, optimal_codons):
    """Frequency of optimal codons (Ikemura 1981).

    ``optimal_codons`` is the set of codons judged optimal for the host, usually
    those with RSCU > 1 in a highly expressed reference set.
    """
    opt = tot = 0
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa is None or aa == "*" or len(FAMILY[aa]) == 1:
            continue
        tot += 1
        if c in optimal_codons:
            opt += 1
    return opt / tot if tot else None


def scuo(seq):
    """Synonymous codon usage order, an entropy measure (Wan et al. 2004).

    0 means codons within each family are used evenly; 1 means a single codon is
    used throughout. Unlike ENC it is bounded and behaves sensibly on short
    sequences, which makes it the safer index for short resistance genes.
    """
    counts = defaultdict(lambda: defaultdict(int))
    for c in codon_list(seq):
        aa = CODON2AA.get(c)
        if aa and aa != "*":
            counts[aa][c] += 1
    num = den = 0.0
    for aa, fam in DEGENERATE.items():
        obs = counts.get(aa)
        if not obs:
            continue
        n = sum(obs.values())
        if n < 1:
            continue
        k = len(fam)
        h = -sum((v / n) * math.log(v / n) for v in obs.values() if v)
        h_max = math.log(k)
        o = (h_max - h) / h_max if h_max else 0.0
        num += o * n
        den += n
    return num / den if den else None
