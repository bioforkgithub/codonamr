"""Synonymous randomisation designs.

These generate null sequences that hold some properties fixed while varying the
one under test. Which properties a null preserves decides what a significant
result means, so each design states exactly what it holds constant.

References
----------
Katz L, Burge CB (2003) Genome Res 13:2042-2051.  (dinucleotide-preserving shuffle)
Workman C, Krogh A (1999) Nucleic Acids Res 27:4816-4822.  (why it is needed)
"""
from collections import defaultdict

from .genetic_code import CODON2AA, FAMILY, codon_list

__all__ = ["randomise_replace", "randomise_shuffle", "randomise_shuffle_dinuc",
           "junction_counts"]


def randomise_replace(seq, preferred, p, rng):
    """Design A: replace each codon with a synonymous one.

    The preferred codon of each family is drawn with marginal probability ``p``
    and the remainder spread evenly. ``p = 1`` gives maximal bias (ENC near 20);
    ``p = 1/k`` gives uniform usage (ENC near 61).

    Changes codon composition, and therefore ENC, CAI, GC3 and the peptide's
    nucleotide context. Holds the peptide fixed. Use it to ask what codon
    *composition* does.
    """
    out = []
    for c in codon_list(seq):
        aa = CODON2AA[c]
        fam = FAMILY[aa]
        k = len(fam)
        if k == 1:
            out.append(c)
            continue
        q = (p - 1.0 / k) / (1.0 - 1.0 / k) if p > 1.0 / k else 0.0
        out.append(preferred[aa] if rng.random() < q else fam[rng.randrange(k)])
    return "".join(out)


def randomise_shuffle(seq, rng):
    """Design B: permute codon positions within each synonymous family.

    Peptide, codon composition, ENC, CAI, GC3 and amino-acid order are all
    invariant by construction. Only codon *order* changes. Use it to ask what
    codon position does, for instance to mRNA secondary structure.

    Caution: this does not preserve dinucleotide composition across codon
    junctions, and folding energy depends on base stacking. A result from this
    design can reflect dinucleotide composition rather than codon order. Use
    :func:`randomise_shuffle_dinuc` to separate the two.
    """
    codons = codon_list(seq)
    pos = defaultdict(list)
    for i, c in enumerate(codons):
        if len(FAMILY[CODON2AA[c]]) > 1:
            pos[CODON2AA[c]].append(i)
    out = list(codons)
    for idx in pos.values():
        vals = [codons[i] for i in idx]
        rng.shuffle(vals)
        for i, v in zip(idx, vals):
            out[i] = v
    return "".join(out)


def junction_counts(codons):
    """Dinucleotides spanning codon boundaries: base 3 of codon i, base 1 of i+1.

    Codon-internal dinucleotides are fixed by codon composition alone, so they
    are invariant under any within-family permutation. These junctions are the
    only part of dinucleotide composition that codon order can change.
    """
    d = defaultdict(int)
    for i in range(len(codons) - 1):
        d[codons[i][2] + codons[i + 1][0]] += 1
    return d


def randomise_shuffle_dinuc(seq, rng, max_iter=None, tol=0):
    """Design B with dinucleotide control.

    Gives the Katz & Burge (2003) DicodonShuffle guarantee by constrained search
    rather than Eulerian construction: start from an unconstrained within-family
    permutation, then accept only synonymous swaps that reduce the L1 distance
    between shuffled and native junction-dinucleotide profiles.

    Returns ``(sequence, residual_L1)``. A residual of 0 means dinucleotide
    composition is matched exactly, so any remaining effect is attributable to
    codon order alone. Report the residual; do not assume it reached zero.
    """
    native = codon_list(seq)
    target = junction_counts(native)
    cur = codon_list(randomise_shuffle(seq, rng))
    n = len(cur)
    if n < 3:
        return "".join(cur), 0

    have = junction_counts(cur)
    dist = sum(abs(target.get(k, 0) - have.get(k, 0)) for k in set(target) | set(have))

    pos = defaultdict(list)
    for i, c in enumerate(cur):
        if len(FAMILY[CODON2AA[c]]) > 1:
            pos[CODON2AA[c]].append(i)
    fams = [v for v in pos.values() if len(v) > 1]
    if not fams or dist <= tol:
        return "".join(cur), dist

    if max_iter is None:
        max_iter = 60 * n

    def junctions_at(seq_c, idx):
        out = []
        if idx > 0:
            out.append((idx - 1, seq_c[idx - 1][2] + seq_c[idx][0]))
        if idx < n - 1:
            out.append((idx, seq_c[idx][2] + seq_c[idx + 1][0]))
        return out

    for _ in range(max_iter):
        if dist <= tol:
            break
        fam = fams[rng.randrange(len(fams))]
        i, j = rng.choice(fam), rng.choice(fam)
        if i == j or cur[i] == cur[j]:
            continue
        affected = dict(junctions_at(cur, i) + junctions_at(cur, j))
        cur[i], cur[j] = cur[j], cur[i]
        after = dict(junctions_at(cur, i) + junctions_at(cur, j))
        delta = 0
        for k, old in affected.items():
            new = after[k]
            if old == new:
                continue
            have[old] -= 1
            delta += (abs(target.get(old, 0) - have[old])
                      - abs(target.get(old, 0) - have[old] - 1))
            have[new] += 1
            delta += (abs(target.get(new, 0) - have[new])
                      - abs(target.get(new, 0) - have[new] + 1))
        if delta > 0:
            cur[i], cur[j] = cur[j], cur[i]
            for k, old in affected.items():
                new = after[k]
                if old == new:
                    continue
                have[new] -= 1
                have[old] += 1
        else:
            dist += delta
    return "".join(cur), dist
