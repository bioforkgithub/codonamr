"""Standard genetic code tables and synonymous families.

Codon order follows the TCAG convention used by CodonW and by most codon-usage
literature, so tables printed by this package line up with published ones.
"""
from collections import defaultdict

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON2AA = dict(zip(CODONS, AAS))
STOPS = frozenset(c for c, a in CODON2AA.items() if a == "*")

_fam = defaultdict(list)
for _c, _a in CODON2AA.items():
    if _a != "*":
        _fam[_a].append(_c)
FAMILY = dict(_fam)

#: the 18 amino acids with more than one codon (excludes Met, Trp and stops)
DEGENERATE = {a: cs for a, cs in FAMILY.items() if len(cs) > 1}

#: amino acids encoded by a single codon, which carry no synonymous information
SINGLETON = {a for a, cs in FAMILY.items() if len(cs) == 1}


def codon_list(seq):
    """Split an in-frame DNA string into codons."""
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
